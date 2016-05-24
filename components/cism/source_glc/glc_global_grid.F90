!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module glc_global_grid

!BOP
! !MODULE: glc_global_grid
!
! !DESCRIPTION:
!  This module contains grid info and routines for setting up the
!  global glc grid quantities.
!  
!  This module is very simple because all we need to do is read in
!   grid size, latitude, and longitude from a file.  This information
!  is needed to receive data from the coupler and pass it to GLINT.  
!
!  Later we may have a more sophisticated, POP-like routine to set up
!  the high-resolution regional ice sheet grid(s).
!
! !REVISION HISTORY:
!  SVN:$Id: grid.F90 808 2006-04-28 17:06:38Z njn01 $
!  Adapted by William Lipscomb from grid.F90 in POP
!
! !USES:

   use glc_kinds_mod
   use glc_communicate
   use glc_broadcast
   use glc_constants
   use glc_exit_mod
   use glint_global_grid
   
   use shr_sys_mod, only : shr_sys_flush
   use shr_file_mod, only : shr_file_getunit, shr_file_freeunit

   implicit none
   private
   save

   public  :: init_glc_grid, glc_landmask, glc_landfrac, glc_grid

! !PUBLIC DATA MEMBERS:

   ! Note that glc_grid, glc_landmask and glc_landfrac only have valid data on the master task
   ! On other tasks, these variables are allocated, but with size 0
   
   type(global_grid) ::   &
      glc_grid        ! info (nx, ny, lat, lon, area) for coupling grid (indexed S to N)

   integer, dimension(:,:), allocatable ::  &
      glc_landmask          ! landmask = 1 where global grid has valid data from CLM
                            ! (i.e., gridcells with nonzero land area)

   real(r8), dimension(:,:), allocatable ::  &
      glc_landfrac          ! land fraction, between 0 and 1
 
!EOP
!BOC
!-----------------------------------------------------------------------
!
!  module private data
!
!-----------------------------------------------------------------------

   character (char_len) ::  &
      horiz_grid_opt       ,&! option for getting horiz grid info
      mask_varname         ,&! variable name for landmask
      frac_varname           ! variable name for landfrac

   character (char_len_long) ::  &
      horiz_grid_file        ! input file for reading horiz grid info

!EOC
!***********************************************************************

 contains

!***********************************************************************

!BOP
! !IROUTINE: init_glc_grid
! !INTERFACE:

 subroutine init_glc_grid

! !DESCRIPTION:
!  Initialize grid quantities
!
! !USES:
   use glc_files, only : nml_filename
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   namelist /grid_nml/ horiz_grid_opt, horiz_grid_file, mask_varname, frac_varname

   integer (int_kind) :: &
      nml_error           ! namelist i/o error flag

!-----------------------------------------------------------------------
!
!  read input namelist for grid setup options
!
!-----------------------------------------------------------------------

   horiz_grid_opt       = 'unknown_horiz_grid_opt'
   horiz_grid_file      = 'unknown_horiz_grid_file'
   mask_varname         = 'unknown mask varname'
   frac_varname         = 'unknown frac varname'

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      do while (nml_error > 0)
         read(nml_in, nml=grid_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
      call exit_glc(sigAbort,'ERROR reading grid_nml')
   endif

   if (my_task == master_task) then
      write(stdout,blank_fmt)
      write(stdout,ndelim_fmt)
      write(stdout,blank_fmt)
      write(stdout,*) ' Grid:'
      write(stdout,blank_fmt)
      write(stdout,*) ' grid_nml namelist settings:'
      write(stdout,blank_fmt)
      write(stdout, grid_nml)
   endif

   call broadcast_scalar(horiz_grid_opt,     master_task)
   call broadcast_scalar(horiz_grid_file,    master_task)
   call broadcast_scalar(mask_varname,       master_task)
   call broadcast_scalar(frac_varname,       master_task)

   write(*,*) my_task, frac_varname
!-----------------------------------------------------------------------
!
!  output grid setup options to log file
!
!-----------------------------------------------------------------------

   if (my_task == master_task) then
      write(stdout,delim_fmt)
      write(stdout,blank_fmt)
      write(stdout,'(a13)') ' Grid options'
      write(stdout,blank_fmt)
   endif

!-----------------------------------------------------------------------
!
!  set up horizontal grid
!
!-----------------------------------------------------------------------

   select case (horiz_grid_opt)
   case ('internal')
      if (my_task == master_task) then
         write(stdout,'(a36)') ' Creating horizontal grid internally'
      endif
      call horiz_grid_internal
   case ('file')
      if (my_task == master_task) then
         write(stdout,*) 'Reading horizontal grid from file:', &
                          trim(horiz_grid_file)
         call shr_sys_flush(stdout)
      endif
      call read_horiz_grid(horiz_grid_file, mask_varname, frac_varname)
   case default
      call exit_glc(sigAbort,'ERROR: unknown horizontal grid option')
   end select

!-----------------------------------------------------------------------
!EOC

 call flushm (stdout)

 end subroutine init_glc_grid

!***********************************************************************
!BOP
! !IROUTINE: horiz_grid_internal
! !INTERFACE:

 subroutine horiz_grid_internal

! !DESCRIPTION:
!  Creates a lat/lon grid with equal spacing in each direction
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i,j,ig,jg,jm1,n    ! dummy counters

   real (r8) :: &
      dlat, dlon,       &! lat/lon spacing for idealized grid
      lathalf,          &! lat at T points
      xdeg               ! temporary longitude variable

!lipscomb - TO DO - Enable internal grid generation, if desired
   call exit_glc(sigAbort, 'Internal grid generation not yet enabled')

!-----------------------------------------------------------------------
!EOC

 end subroutine horiz_grid_internal

!***********************************************************************
!BOP
! !IROUTINE: read_horiz_grid
! !INTERFACE:

 subroutine read_horiz_grid(horiz_grid_file, mask_varname, frac_varname)

! !DESCRIPTION:
!  Reads horizontal grid, landmask, and landfrac from input grid file
!  Note that the data are only valid on the master task - other tasks still create a
!  glc_grid, glc_landmask and glc_landfrac, but with size 0
!
! !REVISION HISTORY:
!  same as module

! !USES:

   use glint_example_clim, only: read_ncdf, read_ncdf_ggrid
   use glint_global_grid

! !INPUT PARAMETERS:

   character (*), intent(in) :: &
      horiz_grid_file  ,&! filename of file containing grid data
      mask_varname     ,&! variable name for landmask
      frac_varname       ! variable name for landfrac
!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer(i4) :: i, j            ! indices
   integer(i4) :: nx, ny          ! global grid dimensions

   real(r8), dimension(:,:), pointer ::  &
      landmask => null()      ,&! landmask, glc/clm grid
      landfrac => null()        ! landfrac, glc/clm grid

   real(r8) :: &
      latn, lats, lone, lonw   ! lat and lon of cell edges (radians)

!-----------------------------------------------------------------------
!
! Extract global grid and mask information from netCDF topography file.
! For now we use GLINT subroutines to do this.
!
!-----------------------------------------------------------------------

! to do - Make sure the read is done correctly
!!   call read_ncdf(horiz_grid_file, mask_varname, landmask, glc_grid)
!!   call read_ncdf(horiz_grid_file, frac_varname, landfrac, glc_grid)

   ! Initialize the grid and masks

   ! Master task gets full grid, others get 0-size grid
   if (my_task==master_task) then
      call read_ncdf_ggrid(horiz_grid_file, glc_grid)

      call read_ncdf(horiz_grid_file, mask_varname, landmask)
      call read_ncdf(horiz_grid_file, frac_varname, landfrac)
  
   else
      call create_empty_grid
   end if

   nx = glc_grid%nx
   ny = glc_grid%ny

!lipscomb - GLINT assumes the grid is indexed N to S and automatically sets
!            lat_bound(1) = 90, lat_bound(ny+1) = -90.
!           Reverse that convention here.

   if (ny > 0) then
      glc_grid%lat_bound(1)    = -90._r8
      glc_grid%lat_bound(ny+1) =  90._r8
   end if

   do i = 1, nx
      if (glc_grid%lon_bound(i) < 0._r8)   &
           glc_grid%lon_bound(i) = glc_grid%lon_bound(i) + 360._r8
   enddo
 
   ! Set glc_landmask and glc_landfrac

   allocate(glc_landmask(nx,ny))
   allocate(glc_landfrac(nx,ny))

   if (my_task==master_task) then
      glc_landmask(:,:) = nint(landmask(:,:))

      glc_landfrac(:,:) = landfrac(:,:)
   endif

   ! Make sure glc_landmask and glc_landfrac have expected values.

   do j = 1, ny
   do i = 1, nx
 
      if (glc_landmask(i,j) /= 0 .and. glc_landmask(i,j) /= 1) then
         write(stdout,*) 'landmask has invalid value: i, j, landmask =', i, j, glc_landmask(i,j)
         call exit_glc(sigAbort, 'invalid landmask value')
      endif

      if (glc_landfrac(i,j) < 0._r8 .or. glc_landfrac(i,j) > 1._r8) then
         write(stdout,*) 'landfrac has invalid value: i, j, landfrac =', i, j, glc_landfrac(i,j)
         call exit_glc(sigAbort, 'invalid landfrac value')
      endif

   enddo
   enddo

   ! compute grid cell area
   ! Note: Global grid is indexed from south to north, so the south edge of cell (i,j+1)
   !       is the north edge of cell (i,j)

   allocate(glc_grid%box_areas(nx,ny))

   do j = 1, ny
   do i = 1, nx

      latn = glc_grid%lat_bound(j+1) * pi/180._r8   ! degrees to radians
      lats = glc_grid%lat_bound(j)   * pi/180._r8
      latn = pi/2._r8 - latn  ! so lat = 0 at NP, = pi at SP
      lats = pi/2._r8 - lats  
      lone = glc_grid%lon_bound(i+1) * pi/180._r8
      lonw = glc_grid%lon_bound(i)   * pi/180._r8
      if (lone < lonw) lone = lone + 2._r8*pi
      glc_grid%box_areas(i,j) = radius**2 * (cos(latn)-cos(lats)) * (lone-lonw)

      ! Make sure area is positive
      if (glc_grid%box_areas(i,j) <= 0._r8) then
         if (verbose) then
            write(stdout,*) 'Negative area: i, j, area =', i, j, glc_grid%box_areas(i,j)
            write(stdout,*) 'latn, lats =', latn, lats
            write(stdout,*) 'cos(latn), cos(lats) =', cos(latn), cos(lats)
            write(stdout,*) 'lone, lonw =', lone, lonw
            write(stdout,*) 'latb(j), latb(j+1) =', glc_grid%lat_bound(j), &
                                                    glc_grid%lat_bound(j+1)
            write(stdout,*) 'lonb(i), lonb(i+1) =', glc_grid%lon_bound(i), &
                                                    glc_grid%lon_bound(i+1)
         endif
         call exit_glc(sigAbort, 'Negative gridcell area on glc grid')
      endif

   enddo
   enddo

   if (verbose .and. my_task==master_task) then
      write(stdout,*) ''
      write(stdout,*) 'Horizontal grid:'
      write(stdout,*) 'nx =', nx
      write(stdout,*) 'ny =', ny
!      write(stdout,*) 'lats =', glc_grid%lats(:)  
!      write(stdout,*) 'lons =', glc_grid%lons(:)
!      write(stdout,*) 'lat_bound =', glc_grid%lat_bound(:)
!      write(stdout,*) 'lon_bound =', glc_grid%lon_bound(:)
      write(stdout,*) ''
      !TODO - Make sure iglint_global and jglint_global are defined appropriately for the global grid
      !      (Currently hardwired in glint_type)
      ! i = iglint_global
      ! j = jglint_global   ! N to S global indexing as in Glint
      ! write(stdout,*) 'Test point, i, j, =', i, j
      ! write(stdout,*) 'lat, lon =', glc_grid%lats(j), glc_grid%lons(i)
      ! write(stdout,*) 'area =', glc_grid%box_areas(i,j)
      ! write(stdout,*) 'landmask =', glc_landmask(i,j)
      ! write(stdout,*) 'landfrac =', glc_landfrac(i,j)
      ! write(stdout,*) 'frac of earth =', glc_grid%box_areas(i,j) / (4._r8*pi*radius*radius)
      write(stdout,*) 'Leaving read_horiz_grid'
      call shr_sys_flush(stdout)
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine read_horiz_grid


!***********************************************************************
!BOP
! !IROUTINE: create_empty_grid
! !INTERFACE:

 subroutine create_empty_grid

! !DESCRIPTION:
!  Sets up glc_grid as a 0-size grid
!
! !REVISION HISTORY:
!  Created by Bill Sacks, Jan 18, 2013

! !USES:

   use glint_global_grid, only : new_global_grid

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   real(r8), dimension(:), allocatable :: lons, lats

!-----------------------------------------------------------------------
!
!  begin code
!
!-----------------------------------------------------------------------

   allocate(lons(0))
   allocate(lats(0))

   call new_global_grid(glc_grid, lons, lats)

   deallocate(lons)
   deallocate(lats)

!-----------------------------------------------------------------------
!EOC

 end subroutine create_empty_grid

!***********************************************************************

 end module glc_global_grid

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

