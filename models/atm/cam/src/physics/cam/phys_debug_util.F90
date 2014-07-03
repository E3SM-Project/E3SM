module phys_debug_util

!----------------------------------------------------------------------------------------

! Module to facilitate debugging of physics parameterizations.
!
! The user requests a location for debugging in lat/lon coordinates
! (degrees).  The initialization routine does a global search to find the
! column in the physics grid closest to the requested location.  The local
! indices of that column in the physics decomposition are stored as module
! data.  The user code then passes the local chunk index of the chunked
! data into the subroutine that will write diagnostic information for the
! column.  The function phys_debug_col returns the local column index if
! the column of interest is contained in the chunk, and zero otherwise.
! Printing is done only if a column index >0 is returned.
!
! Phil Rasch, B. Eaton, Feb 2008
!----------------------------------------------------------------------------------------

use shr_kind_mod,  only: r8 => shr_kind_r8
use phys_grid,     only: phys_grid_find_col, get_rlat_p, get_rlon_p
use spmd_utils,    only: masterproc, iam
use cam_logfile,   only: iulog
use abortutils,    only: endrun

implicit none
private
save

real(r8), parameter :: uninit_r8 = huge(1._r8)

! Public methods
public phys_debug_readnl  ! read namelist input
public phys_debug_init    ! initialize the method to a chunk and column
public phys_debug_col     ! return local column index in debug chunk

! Namelist variables
real(r8) :: phys_debug_lat = uninit_r8 ! latitude of requested debug column location in degrees
real(r8) :: phys_debug_lon = uninit_r8 ! longitude of requested debug column location in degrees


integer :: debchunk = -999            ! local index of the chuck we will debug
integer :: debcol   = -999            ! the column within the chunk we will debug

!================================================================================
contains
!================================================================================

subroutine phys_debug_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'phys_debug_readnl'

   namelist /phys_debug_nl/ phys_debug_lat, phys_debug_lon
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'phys_debug_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, phys_debug_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(phys_debug_lat, 1, mpir8, 0, mpicom)
   call mpibcast(phys_debug_lon, 1, mpir8, 0, mpicom)
#endif

end subroutine phys_debug_readnl

!================================================================================

subroutine phys_debug_init()

   integer  :: owner, lchunk, icol
   real(r8) :: deblat, deblon
   !-----------------------------------------------------------------------------

   ! If no debug column specified then do nothing
   if (phys_debug_lat == uninit_r8 .or. phys_debug_lon == uninit_r8) return

   ! User has specified a column location for debugging.  Find the closest
   ! column in the physics grid.
   call phys_grid_find_col(phys_debug_lat, phys_debug_lon, owner, lchunk, icol)

   ! If the column is owned by this process then save its local indices
   if (iam == owner) then
      debchunk         = lchunk
      debcol           = icol
      deblat           = get_rlat_p(lchunk, icol)*57.296_r8  ! approximate conversion for log output only
      deblon           = get_rlon_p(lchunk, icol)*57.296_r8
      write(iulog,*) 'phys_debug_init: debugging column at lat=', deblat, '  lon=', deblon
   end if

end subroutine phys_debug_init

!================================================================================

integer function phys_debug_col(chunk)

   integer,  intent(in) :: chunk
   !-----------------------------------------------------------------------------

   if (chunk == debchunk) then
      phys_debug_col = debcol
   else
      phys_debug_col = 0
   endif

end function phys_debug_col

!================================================================================

end module phys_debug_util
