!-----------------------------------------------------------------------
module massfix
!----------------------------------------------------------------------- 
! 
! Purpose: Module for mass fixer, contains global integrals
!
!-----------------------------------------------------------------------
!
! Written by: Dani Bundy Coleman, Oct 2004 
! 
!
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use constituents, only: pcnst

!-----------------------------------------------------------------------
   implicit none
!
! By default everything is private to this module
!
   private
!
! Public interfaces
!

   public hw1, hw2, hw3, alpha         ! Needs to be public for restart

!
! Module data
!
   real(r8) :: hw1(pcnst)              ! Pre-SLT global integral of constituent
   real(r8) :: hw2(pcnst)              ! Post-SLT global integral of const.
   real(r8) :: hw3(pcnst)              ! Global integral for denom. of expr. for alpha
   real(r8) :: alpha(pcnst)            ! alpha(m) = ( hw1(m) - hw2(m) )/hw3(m)


end module massfix
