!=======================================================================
!
!BOP
!
! !MODULE: ice_work - globally accessible, temporary work arrays
!
! !DESCRIPTION:
!
! Declare globally accessible, temporary work arrays to conserve memory.
! These arrays should be used only within a single subroutine!
!
! !REVISION HISTORY:
!  SVN:$Id: ice_work.F90 37 2006-11-29 18:06:44Z eclare $
!
! authors Elizabeth C. Hunke and William H. Lipscomb, LANL
!
! 2004: Block structure added by William Lipscomb
! 2006: Converted to free source form (F90) by Elizabeth Hunke
!
! !INTERFACE:
!
      module ice_work
!
! !USES:
!
      use ice_kinds_mod
      use ice_blocks
      use ice_domain_size
!
!EOP
!
      implicit none

      ! global

      real (kind=dbl_kind), dimension(:,:), allocatable :: &
         work_g1, &
         work_g2, &
         work_g3

      real (kind=real_kind), dimension(:,:), allocatable :: &
         work_gr

      real (kind=real_kind), dimension(:,:,:), allocatable :: &
         work_gr3

      integer(kind=int_kind), dimension(:,:), allocatable :: &
         work_gi4

      integer(selected_int_kind(13)), dimension(:,:), allocatable :: &
         work_gi8

      ! all blocks
      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         work1, &
         work2

      ! local (single block)
      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         worka, &
         workb, &
         workc, &
         workd

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !ROUTINE: init_work - initialize work arrays
!
! !DESCRIPTION:
!
! Initialize work arrays
!
! !REVISION HISTORY:
!
! author: William H. Lipscomb, LANL
!
! !INTERFACE:
!
      subroutine init_work
!
! !USES:
!
      use ice_constants
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      work1(:,:,:) = c0
      work2(:,:,:) = c0

      worka(:,:) = c0
      workb(:,:) = c0
      workc(:,:) = c0
      workd(:,:) = c0

      end subroutine init_work

!=======================================================================

      end module ice_work

!=======================================================================
