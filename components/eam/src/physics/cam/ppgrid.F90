
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
  public nvar_dirOA
  public nvar_dirOL
  public indexb

! Grid point resolution parameters

#ifdef PPCOLS
   integer pcols      ! max number of columns per chunk (set at compile-time)
#endif
   integer psubcols   ! number of sub-columns (max)
   integer pver       ! number of vertical levels
   integer pverp      ! pver + 1
   !added for ogwd
   integer nvar_dirOA
   integer nvar_dirOL
   integer indexb

#ifdef PPCOLS
   parameter (pcols     = PCOLS)
#endif
   parameter (psubcols  = PSUBCOLS)
   parameter (pver      = PLEV)
   parameter (pverp     = pver + 1  )
   !added for ogwd
   parameter (nvar_dirOA =2+1 )!avoid bug when nvar_dirOA is 2
   parameter (nvar_dirOL =180)!set for 360 degrees wind direction
   parameter (indexb = 3232)!set for 3km-inputs
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
