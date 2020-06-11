
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


! Grid point resolution parameters

#ifdef PPCOLS
   integer pcols      ! max number of columns per chunk (set at compile-time)
#endif
   integer psubcols   ! number of sub-columns (max)
   integer pver       ! number of vertical levels
   integer pverp      ! pver + 1

#ifdef PPCOLS
   parameter (pcols     = PCOLS)
#endif
   parameter (psubcols  = PSUBCOLS)
   parameter (pver      = PLEV)
   parameter (pverp     = pver + 1  )
!
! start, end indices for chunks owned by a given MPI task
! (set in phys_grid_init).
!
   integer :: begchunk = 0            ! 
   integer :: endchunk = -1           ! 

#ifndef PPCOLS
!
! pcols set in physgrid_init
   integer :: pcols = -1    ! max number of columns per chunk (set at runtime)
#endif

end module ppgrid
