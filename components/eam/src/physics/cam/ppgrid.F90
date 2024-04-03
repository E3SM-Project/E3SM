
module ppgrid

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Initialize physics grid resolution parameters
!  for a chunked data structure
! 
! Author: 
! 
!-----------------------------------------------------------------------

  implicit none
  private
  save
 
  public begchunk
  public endchunk
  public pcols
  public psubcols
  public pver
  public pverp
!!====Jinbo Xie====
  public nvar_dirOA
  public nvar_dirOL
  public indexb
!!====Jinbo Xie====


! Grid point resolution parameters
!!Jinbo Xie
#ifdef PPCOLS
   integer pcols      ! max number of columns per chunk (set at compile-time)
#endif
   !!Jinbo Xie
   integer psubcols   ! number of sub-columns (max)
   integer pver       ! number of vertical levels
   integer pverp      ! pver + 1
   !====Jinbo Xie====
   integer nvar_dirOA
   integer nvar_dirOL
   integer indexb
   !====Jinbo Xie====
   !!Jinbo Xie
#ifdef PPCOLS
   parameter (pcols     = PCOLS)
#endif
   !!Jinbo Xie
   parameter (psubcols  = PSUBCOLS)
   parameter (pver      = PLEV)
   parameter (pverp     = pver + 1  )
   !!====Jinbo Xie====
   parameter (nvar_dirOA =2+1 )!Jinbo Xie avoid bug when nvar_dirOA is 2
   parameter (nvar_dirOL =180)
   parameter (indexb = 3232)
   !!====Jinbo Xie====
!
! start, end indices for chunks owned by a given MPI task
! (set in phys_grid_init).
!
   integer :: begchunk = 0            ! 
   integer :: endchunk = -1           ! 


   !!Jinbo Xie
#ifndef PPCOLS
! pcols set in physgrid_init
   integer :: pcols = -1    ! max number of columns per chunk (set at runtime)
#endif
   !!Jinbo Xie

end module ppgrid
