      module parkind

      implicit none
      save

!------------------------------------------------------------------
! rrtmg kinds
! Define integer and real kinds for various types.
!
! Initial version: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!------------------------------------------------------------------

!
!     integer kinds
!     -------------
!
      integer, parameter :: kind_ib = selected_int_kind(13)  ! 8 byte integer
      integer, parameter :: kind_im = selected_int_kind(6)   ! 4 byte integer
      integer, parameter :: kind_in = kind(1)                ! native integer

!
!     real kinds
!     ----------
!
      integer, parameter :: kind_rb = selected_real_kind(12) ! 8 byte real
      integer, parameter :: kind_rm = selected_real_kind(6)  ! 4 byte real
      integer, parameter :: kind_rn = kind(1.0)              ! native real

      end module parkind
