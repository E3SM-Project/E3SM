      module rrlw_tbl

      use parkind, only : im => kind_im, rb => kind_rb

      implicit none
      save

!------------------------------------------------------------------
! rrtmg_lw exponential lookup table arrays

! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, Jun 2006
! Revised: MJIacono, AER, Aug 2007
! Revised: MJIacono, AER, Aug 2008
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

      integer(kind=im), parameter :: ntbl = 10000

      real(kind=rb), parameter :: tblint = 10000.0_rb

      real(kind=rb) , dimension(0:ntbl) :: tau_tbl
      real(kind=rb) , dimension(0:ntbl) :: exp_tbl
      real(kind=rb) , dimension(0:ntbl) :: tfn_tbl

      real(kind=rb), parameter :: pade = 0.278_rb
      real(kind=rb) :: bpade

      end module rrlw_tbl

