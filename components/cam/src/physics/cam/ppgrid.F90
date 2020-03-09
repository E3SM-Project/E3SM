
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
  public ppcols
  public psubcols
  public pver
  public pverp


! Grid point resolution parameters

   integer ppcols     ! number of columns (default for runtime-set max pcols)
   integer psubcols   ! number of sub-columns (max)
   integer pver       ! number of vertical levels
   integer pverp      ! pver + 1

   parameter (ppcols    = PCOLS)
   parameter (psubcols  = PSUBCOLS)
   parameter (pver      = PLEV)
   parameter (pverp     = pver + 1  )
!
! start, end indices for chunks owned by a given MPI task
! (set in phys_grid_init).
!
   integer :: begchunk = 0            ! 
   integer :: endchunk = -1           ! 

!
! pcols set in physgrid_init
   integer :: pcols = PCOLS    ! number of columns (runtime-set max)

end module ppgrid
