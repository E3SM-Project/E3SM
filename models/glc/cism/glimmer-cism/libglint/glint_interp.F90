!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glint_interp.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2013
!   Glimmer-CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of Glimmer-CISM.
!
!   Glimmer-CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   Glimmer-CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with Glimmer-CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

module glint_interp

  !*FD Downscaling and upscaling routines for use in Glint

  use glimmer_global, only: dp, sp
  use glimmer_map_types
  use glint_mpinterp
  use glimmer_paramets, only: stdout, GLC_DEBUG

  implicit none

  type downscale

     !*FD Derived type containing indexing 
     !*FD information for downscaling. This type was 
     !*FD included for speed. Four of the arrays contained in it
     !*FD are arrays of the indices of the corners 
     !*FD of the global grid-boxes within which the given
     !*FD local grid point lies.

     real(dp),dimension(:,:),pointer :: llats => null() !*FD The latitude of each point in x-y space.
     real(dp),dimension(:,:),pointer :: llons => null() !*FD The longitude of each point in x-y space.

     integer, dimension(:,:,:),pointer :: xloc => null() !*FD The x-locations of the corner points of the
                                                         !*FD interpolation domain.
     integer, dimension(:,:,:),pointer :: yloc => null() !*FD The y-locations of the corner points of the
                                                         !*FD interpolation domain.
     integer, dimension(:,:), pointer :: xin => null() !*FD x-locations of global cell the point is in
     integer, dimension(:,:), pointer :: yin => null() !*FD y-locations of global cell the point is in

     real(dp),dimension(:,:),  pointer :: xfrac => null()
     real(dp),dimension(:,:),  pointer :: yfrac => null()
     real(dp),dimension(:,:),pointer :: sintheta => NULL()  !*FD sines of grid angle relative to north.
     real(dp),dimension(:,:),pointer :: costheta => NULL()  !*FD coses of grid angle relative to north.
     type(mpinterp) :: mpint !*FD Parameters for mean-preserving interpolation
     logical :: use_mpint = .false. !*FD set true if we're using mean-preserving interpolation
     integer,dimension(:,:),pointer :: lmask => null()  !*FD mask = 1 where downscaling is valid 
                                                        !*FD mask = 0 elsewhere

  end type downscale

  type upscale

     !*FD Derived type containing indexing information
     !*FD for upscaling by areal averaging.

     integer, dimension(:,:),pointer :: gboxx => null() !*FD $x$-indices of global grid-box 
                                                        !*FD containing given local grid-box.
     integer, dimension(:,:),pointer :: gboxy => null() !*FD $y$-indices of global grid-box 
                                                        !*FD containing given local grid-box.
     integer, dimension(:,:),pointer :: gboxn => null() !*FD Number of local grid-boxes 
                                                        !*FD contained in each global box.
     logical                         :: set = .false.   !*FD Set if the type has been initialised.

  end type upscale

contains

  subroutine new_downscale(downs,proj,ggrid,lgrid,mpint)

    use glint_global_grid
    use glimmer_map_trans
    use glimmer_map_types
    use glimmer_coordinates

    !*FD Initialises a downscale variable,
    !*FD according to given projected and global grids

    ! Arguments

    type(downscale),intent(out)      :: downs   !*FD Downscaling variable to be set
    type(glimmap_proj),intent(in)    :: proj    !*FD Projection to use
    type(global_grid),intent(in)     :: ggrid   !*FD Global grid to use
    type(coordsystem_type),intent(in) :: lgrid  !*FD Local (ice) grid 
    logical,optional :: mpint !*FD Set true if we're using mean-preserving interp

    ! Internal variables

    real(dp) :: llat,llon
    integer :: i,j
    type(upscale) :: ups
    integer,dimension(:,:),pointer :: upsm

    upsm => null()
    ! Allocate arrays

    allocate(downs%xloc(lgrid%size%pt(1),lgrid%size%pt(2),4))
    allocate(downs%yloc(lgrid%size%pt(1),lgrid%size%pt(2),4))
    call coordsystem_allocate(lgrid,downs%xfrac)
    call coordsystem_allocate(lgrid,downs%yfrac)
    call coordsystem_allocate(lgrid,downs%llons)
    call coordsystem_allocate(lgrid,downs%llats)
    call coordsystem_allocate(lgrid,downs%sintheta)
    call coordsystem_allocate(lgrid,downs%costheta)
    call coordsystem_allocate(lgrid,downs%xin)
    call coordsystem_allocate(lgrid,downs%yin)
    call coordsystem_allocate(lgrid,upsm)
    call coordsystem_allocate(lgrid,downs%lmask)

    ! index local boxes

    call index_local_boxes(downs%xloc,  downs%yloc,   &
                           downs%xfrac, downs%yfrac,  &
                           ggrid, proj, lgrid,        &
                           downs%lmask )

    ! Calculate grid angle

    call calc_grid_angle(downs,proj,lgrid)

    ! Find lats and lons

    do i=1,lgrid%size%pt(1)
       do j=1,lgrid%size%pt(2)
          call glimmap_xy_to_ll(llon,llat,real(i,kind=dp),real(j,kind=dp),proj,lgrid)
          downs%llons(i,j)=llon
          downs%llats(i,j)=llat
       end do
    end do

    ! Initialise mean-preserving interpolation if necessary
    if (present(mpint)) then
       if (mpint) then
          call new_mpinterp(downs%mpint,ggrid)
          downs%use_mpint = .true.
       end if
    end if

    call new_upscale(ups,ggrid,proj,upsm,lgrid)
    downs%xin = ups%gboxx
    downs%yin = ups%gboxy
    deallocate(upsm)

  end subroutine new_downscale

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine interp_wind_to_local(lgrid_fulldomain,zonwind,merwind,downs,xwind,ywind)

    !*FD Interpolates a global wind field 
    !*FD (or any vector field) onto a given projected grid.
    !
    !*FD Currently doesn't work with multiple tasks

    use glimmer_utils
    use glimmer_coordinates
    use glimmer_log
    use parallel, only : tasks

    ! Argument declarations

    type(coordsystem_type), intent(in)     :: lgrid_fulldomain !*FD Target grid on the full domain (i.e., across all tasks)
    real(dp),dimension(:,:),intent(in)     :: zonwind          !*FD Zonal component (input)
    real(dp),dimension(:,:),intent(in)     :: merwind          !*FD Meridional components (input)
    type(downscale),        intent(inout)  :: downs            !*FD Downscaling parameters
    real(dp),dimension(:,:),intent(out)    :: xwind,ywind      !*FD x and y components on the projected grid (output)

    ! Declare two temporary arrays to hold the interpolated zonal and meridional winds

    real(dp),dimension(size(xwind,1),size(xwind,2)) :: tempzw,tempmw

    ! Check input arrays are conformal to one another

    call check_conformal(zonwind,merwind,'interp_wind 1')
    call check_conformal(xwind,ywind,'interp_wind 2')

    ! Interpolate onto the projected grid

    call interp_to_local(lgrid_fulldomain,zonwind,downs,localdp=tempzw)
    call interp_to_local(lgrid_fulldomain,merwind,downs,localdp=tempmw)

    ! Apply rotation

    ! WJS (1-15-13): The following code won't work currently if there is more than 1 task,
    ! because the downs variable applies to the full (non-decomposed) domain, and is only
    ! valid on the master task
    if (tasks > 1) then
       call write_log('interp_wind_to_local only works with a single task', &
                      GM_FATAL, __FILE__, __LINE__)
    end if
    xwind=tempzw*downs%costheta-tempmw*downs%sintheta
    ywind=tempzw*downs%sintheta+tempmw*downs%costheta

  end subroutine interp_wind_to_local

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine interp_to_local (lgrid_fulldomain, global,      downs,   &
                              localsp,          localdp,              &
                              global_fn,        z_constrain,          &
                              gmask,            maskval)

    !*FD Interpolate a global scalar field onto a projected grid. 
    !*FD 
    !*FD This uses a simple bilinear interpolation, which assumes
    !*FD that the global grid boxes are rectangular - i.e. it works
    !*FD in lat-lon space.
    !*FD
    !*FD Either localsp or localdp must be present (or both), depending
    !*FD which precision output is required.
    !*FD
    !*FD Variables referring to the global domain (global, downs,
    !*FD gmask) only need to be valid on the main task

! Cell indexing for (xloc,yloc) is as follows:
!
!       4---------3
!       |         |
!       |         |
!       |         |
!       1---------2
!

    use glimmer_utils
    use glimmer_coordinates
    use glimmer_log
    use parallel, only : main_task, distributed_scatter_var, parallel_halo

    ! Argument declarations

    type(coordsystem_type),  intent(in)           :: lgrid_fulldomain !*FD Local grid, spanning the full domain (across all tasks)
    real(dp), dimension(:,:),intent(in)           :: global           !*FD Global field (input)
    type(downscale),         intent(inout)        :: downs            !*FD Downscaling parameters
    !TODO - Remove localsp as optional argument?
    real(sp),dimension(:,:), intent(out),optional :: localsp          !*FD Local field on projected grid (output) sp
    real(dp),dimension(:,:), intent(out),optional :: localdp          !*FD Local field on projected grid (output) dp
    real(dp),optional,external                    :: global_fn        !*FD Function returning values in global field. This  
                                                                      !*FD may be used as an alternative to passing the
                                                                      !*FD whole array in \texttt{global} if, for instance the
                                                                      !*FD data-set is in a large file, being accessed point by point.
                                                                      !*FD In these circumstances, \texttt{global}
                                                                      !*FD may be of any size, and its contents are irrelevant.
    logical,optional :: z_constrain
    integer, dimension(:,:), intent(in),optional  :: gmask            !*FD = 1 where global data are valid, else = 0
    real(dp), intent(in), optional                :: maskval          !*FD Value to write for masked-out cells 

    ! Local variable declarations

    !TODO - Remove localsp_fulldomain?
    real(sp), dimension(:,:), allocatable :: localsp_fulldomain  ! localsp spanning full domain (all tasks)
    real(dp), dimension(:,:), allocatable :: localdp_fulldomain  ! localdp spanning full domain (all tasks)
    integer  :: i,j                          ! Counter variables for main loop
    real(dp),dimension(4) :: f               ! Temporary array holding the four points in the 
                                             ! interpolation domain.
    real(dp), dimension(size(global,1),size(global,2)) :: g_loc
    logical,  dimension(size(global,1),size(global,2)) :: zeros
    logical :: zc

    integer :: x1, x2, x3, x4
    integer :: y1, y2, y3, y4

    if (present(z_constrain)) then
       zc=z_constrain
    else
       zc=.false.
    end if

    ! check we have one output at least...

    if (.not. (present(localsp) .or. present(localdp)) ) then
       call write_log('Interp_to_local has no output',GM_WARNING,__FILE__,__LINE__)
    endif

    ! Allocate variables to hold result of interpolation
    ! We allocate size 0 arrays on non-main task (rather than leaving variables
    ! unallocated there), because distributed_scatter_var tries to do a deallocate on all tasks
    ! Note that coordsystem_allocate can't be used here because it only works on pointer
    ! variables, and the *_fulldomain variables are non-pointers (as is required for distributed_scatter_var)

    if (present(localsp)) then
       if (main_task) then
          allocate(localsp_fulldomain(lgrid_fulldomain%size%pt(1), lgrid_fulldomain%size%pt(2)))
       else
          allocate(localsp_fulldomain(0,0))
       end if
    end if

    if (present(localdp)) then
       if (main_task) then
          allocate(localdp_fulldomain(lgrid_fulldomain%size%pt(1), lgrid_fulldomain%size%pt(2)))
       else
          allocate(localdp_fulldomain(0,0))
       end if
    end if

    ! Do main interpolation work, just on main task

    if (main_task) then

       ! Do stuff for mean-preserving interpolation

       if (downs%use_mpint) then
          call mean_preserve_interp(downs%mpint,global,g_loc,zeros)
       end if

       ! Main interpolation loop

       do i=1,lgrid_fulldomain%size%pt(1)
          do j=1,lgrid_fulldomain%size%pt(2)

             ! Compile the temporary array f from adjacent points 

             !TODO - This could be handled more efficiently by precomputing arrays that specify
             !       which neighbor gridcell supplies values in each masked-out global gridcell.

             if (present(gmask) .and. present(maskval)) then

                if (downs%lmask(i,j) == 0) then

                   f(1) = maskval
                   f(2) = maskval
                   f(3) = maskval
                   f(4) = maskval

                else

                   x1 = downs%xloc(i,j,1); y1 = downs%yloc(i,j,1)
                   x2 = downs%xloc(i,j,2); y2 = downs%yloc(i,j,2)
                   x3 = downs%xloc(i,j,3); y3 = downs%yloc(i,j,3)
                   x4 = downs%xloc(i,j,4); y4 = downs%yloc(i,j,4)

                   ! if a gridcell is masked out, try to assign a value from a
                   ! neighbor that is not masked out

                   if (gmask(x1,y1) /= 0) then
                      f(1) = global(x1,y1)
                   elseif (gmask(x2,y2) /= 0) then
                      f(1) = global(x2,y2)
                   elseif (gmask(x4,y4) /= 0) then
                      f(1) = global(x4,y4)
                   elseif  (gmask(x3,y3) /= 0) then
                      f(1) = global(x3,y3)
                   else
                      f(1) = maskval
                   endif

                   if (gmask(x2,y2) /= 0) then
                      f(2) = global(x2,y2)
                   elseif (gmask(x1,y1) /= 0) then
                      f(2) = global(x1,y1)
                   elseif (gmask(x3,y3) /= 0) then
                      f(2) = global(x3,y3)
                   elseif  (gmask(x4,y4) /= 0) then
                      f(2) = global(x4,y4)
                   else
                      f(2) = maskval
                   endif

                   if (gmask(x3,y3) /= 0) then
                      f(3) = global(x3,y3)
                   elseif (gmask(x4,y4) /= 0) then
                      f(3) = global(x4,y4)
                   elseif (gmask(x2,y2) /= 0) then
                      f(3) = global(x2,y2)
                   elseif  (gmask(x1,y1) /= 0) then
                      f(3) = global(x1,y1)
                   else
                      f(3) = maskval
                   endif

                   if (gmask(x4,y4) /= 0) then
                      f(4) = global(x4,y4)
                   elseif (gmask(x3,y3) /= 0) then
                      f(4) = global(x3,y3)
                   elseif (gmask(x1,y1) /= 0) then
                      f(4) = global(x1,y1)
                   elseif  (gmask(x2,y2) /= 0) then
                      f(4) = global(x2,y2)
                   else
                      f(4) = maskval
                   endif

                endif    ! lmask = 0

             else        ! gmask and maskval not present

                if (present(global_fn)) then
                   f(1)=global_fn(downs%xloc(i,j,1),downs%yloc(i,j,1))
                   f(2)=global_fn(downs%xloc(i,j,2),downs%yloc(i,j,2))
                   f(3)=global_fn(downs%xloc(i,j,3),downs%yloc(i,j,3))
                   f(4)=global_fn(downs%xloc(i,j,4),downs%yloc(i,j,4))
                else
                   if (downs%use_mpint) then
                      f(1)=g_loc(downs%xloc(i,j,1),downs%yloc(i,j,1))
                      f(2)=g_loc(downs%xloc(i,j,2),downs%yloc(i,j,2))
                      f(3)=g_loc(downs%xloc(i,j,3),downs%yloc(i,j,3))
                      f(4)=g_loc(downs%xloc(i,j,4),downs%yloc(i,j,4))
                   else
                      f(1)=global(downs%xloc(i,j,1),downs%yloc(i,j,1))
                      f(2)=global(downs%xloc(i,j,2),downs%yloc(i,j,2))
                      f(3)=global(downs%xloc(i,j,3),downs%yloc(i,j,3))
                      f(4)=global(downs%xloc(i,j,4),downs%yloc(i,j,4))
                   end if
                end if

             endif  ! gmask and maskval present

             ! Apply the bilinear interpolation

             if (zc.and.zeros(downs%xin(i,j),downs%yin(i,j)).and.downs%use_mpint) then
                if (present(localsp)) localsp_fulldomain(i,j)=0.0_sp
                if (present(localdp)) localdp_fulldomain(i,j)=0.0_dp
             else
                if (present(localsp)) localsp_fulldomain(i,j)=bilinear_interp(downs%xfrac(i,j),downs%yfrac(i,j),f)
                if (present(localdp)) localdp_fulldomain(i,j)=bilinear_interp(downs%xfrac(i,j),downs%yfrac(i,j),f)
             end if

          enddo
       enddo
    end if  ! main_task

    ! Main task scatters interpolated data from the full domain to the task owning each point
    ! Note that distributed_scatter_var doesn't set halo values, so we need to do a halo
    ! update if it's important to have correct values in the halo cells.

    if (present(localsp)) then
       call distributed_scatter_var(localsp, localsp_fulldomain)
       call parallel_halo(localsp)
    end if
    if (present(localdp)) then
       call distributed_scatter_var(localdp, localdp_fulldomain)
       call parallel_halo(localdp)
    end if

    ! We do NOT deallocate the local*_fulldomain variables here, because the
    ! distributed_scatter_var routines do this deallocation

  end subroutine interp_to_local

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine copy_to_local (lgrid_fulldomain, global, downs, local)

    !*FD Do a simple copy of a global scalar field onto a projected grid.
    !*FD 
    !*FD This copies the value from each global cell into all local cells contained
    !*FD within it.
    !*FD
    !*FD Note that, in contrast to interp_to_local, this routine does not support a gmask.
    !*FD
    !*FD Variables referring to the global domain (global, downs) only need to be valid
    !*FD on the main task.

    use glimmer_coordinates
    use parallel, only : main_task, distributed_scatter_var, parallel_halo

    ! Argument declarations

    type(coordsystem_type),  intent(in)  :: lgrid_fulldomain !*FD Local grid, spanning the full domain (across all tasks)
    real(dp), dimension(:,:),intent(in)  :: global           !*FD Global field (input)
    type(downscale),         intent(in)  :: downs            !*FD Downscaling parameters
    real(dp),dimension(:,:), intent(out) :: local            !*FD Local field on projected grid (output)
    
    ! Local variable declarations

    real(dp), dimension(:,:), allocatable :: local_fulldomain  ! local spanning full domain (all tasks)
    integer :: i,j    ! local indices
    integer :: ig,jg  ! global indices

    if (main_task) then
       allocate(local_fulldomain(lgrid_fulldomain%size%pt(1), lgrid_fulldomain%size%pt(2)))
    else
       allocate(local_fulldomain(0,0))
    end if

    ! Do main copying work, just on main task

    if (main_task) then
       do j=1,lgrid_fulldomain%size%pt(2)
          do i=1,lgrid_fulldomain%size%pt(1)
             ig = downs%xin(i,j)
             jg = downs%yin(i,j)
             local_fulldomain(i,j) = global(ig,jg)
          end do
       end do
    end if

    ! Main task scatters interpolated data from the full domain to the task owning each point
    ! Note that distributed_scatter_var doesn't set halo values, so we need to do a halo
    ! update if it's important to have correct values in the halo cells.

    call distributed_scatter_var(local, local_fulldomain)
    call parallel_halo(local)

    ! We do NOT deallocate local_fulldomain here, because the distributed_scatter_var
    ! routine does this deallocation

  end subroutine copy_to_local

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine mean_to_local(proj,lgrid,ggrid,global,localsp,localdp,global_fn)

    !*FD Average a high-resolution global field onto the projected grid
    !*FD This assumes that the global field is sufficiently high-resolution 
    !*FD compared with the local grid - it just averages the points contained 
    !*FD in each local grid-box.
    !*FD
    !*FD This may not work properly with multiple tasks.

    use glimmer_map_types
    use glimmer_map_trans
    use glimmer_coordinates
    use glimmer_utils
    use glimmer_log
    use glint_global_grid

    ! Argument declarations

    type(glimmap_proj),              intent(in)  :: proj      !*FD Target map projection
    type(coordsystem_type),          intent(in)  :: lgrid     !*FD Local grid information
    type(global_grid),               intent(in)  :: ggrid     !*FD Global grid information
    real(dp),dimension(:,:),         intent(in)  :: global    !*FD Global field (input)
    !TODO - Remove localsp?
    real(sp),dimension(:,:),optional,intent(out) :: localsp   !*FD Local field on projected grid (output) sp
    real(dp),dimension(:,:),optional,intent(out) :: localdp   !*FD Local field on projected grid (output) dp
    real(dp),optional, external                  :: global_fn !*FD Function returning values in global field. This  
                                                              !*FD may be used as an alternative to passing the
                                                              !*FD whole array in \texttt{global} if, for instance the
                                                              !*FD data-set is in a large file, being accessed point by point.
                                                              !*FD In these circumstances, \texttt{global}
                                                              !*FD may be of any size, and its contents are irrelevant.

    integer :: i,j,xbox,ybox
    real(dp) :: lat,lon,x,y
    real(dp),dimension(lgrid%size%pt(1),lgrid%size%pt(2)) :: temp_out
    real(dp),dimension(lgrid%size%pt(1),lgrid%size%pt(2)) :: mean_count

    if (.not.present(global_fn)) then 
       if ((lgrid%size%pt(1)/=size(ggrid%lons)).or.(lgrid%size%pt(2)/=size(ggrid%lats))) then
          call write_log('Size mismatch in interp_to_local',GM_FATAL,__FILE__,__LINE__)
       end if
    end if

    ! check we have one output at least...

    if (.not. (present(localsp) .or. present(localdp))) then
       call write_log('mean_to_local has no output',GM_WARNING,__FILE__,__LINE__)
    endif

    ! Zero some things

    mean_count = 0
    temp_out = 0.d0

    ! Loop over all global points

    do i=1,lgrid%size%pt(1)

       lon=ggrid%lons(i)

       do j=1,lgrid%size%pt(2)

          ! Find location in local coordinates

          lat=ggrid%lats(j)  ! (Have already found lat above)
          call glimmap_ll_to_xy(lon,lat,x,y,proj,lgrid)
          xbox=nint(x)
          ybox=nint(y)

          ! Add to appropriate location and update count

          if (xbox >= 1.and.xbox <= lgrid%size%pt(1).and. &
               ybox >= 1.and.ybox <= lgrid%size%pt(2)) then
             if (present(global_fn)) then
                temp_out(xbox,ybox)=temp_out(xbox,ybox)+global_fn(i,j)*ggrid%box_areas(xbox,ybox)
             else
                temp_out(xbox,ybox)=temp_out(xbox,ybox)+global(i,j)*ggrid%box_areas(xbox,ybox)
             end if
             mean_count(xbox,ybox)=mean_count(xbox,ybox)+ggrid%box_areas(xbox,ybox)
          end if

       end do
    end do

    ! Divide by number of contributing points and copy to output

    if (present(localsp)) localsp = temp_out/real(mean_count,dp)
    if (present(localdp)) localdp = temp_out/real(mean_count,dp)

  end subroutine mean_to_local

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine pointwise_to_global(proj,lgrid,local,lons,lats,global)

    !*FD Upscale to global domain by
    !*FD pointwise sampling.
    !*FD
    !*FD Note that this is the mathematically inverse process of the 
    !*FD \texttt{interp\_to\_local} routine.
    !*FD
    !*FD Note that this probably doesn't work right with multiple tasks

    use glimmer_coordinates
    use glimmer_map_trans

    ! Arguments

    type(glimmap_proj),     intent(in)  :: proj      !*FD Projection to use
    type(coordsystem_type), intent(in)  :: lgrid     !*FD Local grid
    real(dp),dimension(:,:),intent(in)  :: local     !*FD Local field (input)
    real(dp),dimension(:,:),intent(out) :: global    !*FD Global field (output)
    real(dp),dimension(:),  intent(in)  :: lats      !*FD Latitudes of grid-points (degrees)
    real(dp),dimension(:),  intent(in)  :: lons      !*FD Longitudes of grid-points (degrees)

    ! Internal variables

    real(dp),dimension(2,2) :: f
    integer :: nxg,nyg,nxl,nyl,i,j,xx,yy
    real(dp) :: x,y
    real(dp),dimension(size(local,1),size(local,2)) :: tempmask

    nxg=size(global,1) ; nyg=size(global,2)
    nxl=size(local,1)  ; nyl=size(local,2)

    do i=1,nxg
       do j=1,nyg
          call glimmap_ll_to_xy(lons(i),lats(j),x,y,proj,lgrid)
          xx=int(x) ; yy=int(y)
          if (nint(x)<=1.or.nint(x)>nxl-1.or.nint(y)<=1.or.nint(y)>nyl-1) then
             global(i,j)=0.0
          else
             f=local(xx:xx+1,yy:yy+1)*tempmask(xx:xx+1,yy:yy+1)
             global(i,j)=bilinear_interp((x-real(xx))/real(1.0,dp),(y-real(yy))/real(1.0,dp),f)
          endif
       enddo
    enddo

  end subroutine pointwise_to_global

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !WHL - Changed name from mean_to_global to local_to_global_avg

  subroutine local_to_global_avg(ups,local,global,mask)

    !*FD Upscale to global domain by areal averaging.
    !*FD
    !*FD Note that:
    !*FD \begin{itemize}
    !*FD \item \texttt{global} output is only valid on the main task
    !*FD \item \texttt{ups} input only needs to be valid on the main task
    !*FD \item \texttt{gboxx} and \texttt{gboxy} are the same size as \texttt{local_fulldomain}
    !*FD \item \texttt{gboxn} is the same size as \texttt{global}
    !*FD \item This method is \emph{not} the mathematical inverse of the
    !*FD \texttt{interp\_to\_local} routine.
    !*FD \end{itemize}

    use parallel, only : main_task, distributed_gather_var
    use nan_mod , only : NaN

    ! Arguments

    type(upscale),          intent(in)  :: ups    !*FD Upscaling indexing data.
    real(dp),dimension(:,:),intent(in)  :: local  !*FD Data on projected grid (input).
    real(dp),dimension(:,:),intent(out) :: global !*FD Data on global grid (output).
    integer, dimension(:,:),intent(in),optional :: mask !*FD Mask for upscaling

    ! Internal variables

    integer :: nxl_full,nyl_full,i,j
    real(dp),dimension(size(local,1),size(local,2)) :: tempmask
    
    ! values of 'local' and 'tempmask' spanning full domain (all tasks)
    real(dp),dimension(:,:), allocatable            :: local_fulldomain
    real(dp),dimension(:,:), allocatable            :: tempmask_fulldomain
    real(dp),dimension(:,:) ,allocatable            :: ncells

    ! Beginning of code
    
    allocate(ncells(size(global,1), size(global,2)))
    
    global(:,:) = NaN

    if (present(mask)) then
       tempmask = mask
    else
       tempmask = 1.d0
    endif

    ! Gather 'local' and 'tempmask' onto main task, which is the only one that does the regridding

    call distributed_gather_var(local, local_fulldomain)
    call distributed_gather_var(tempmask, tempmask_fulldomain)

    ! Main task does regridding

    if (main_task) then

       nxl_full = size(local_fulldomain,1)
       nyl_full = size(local_fulldomain,2)
       global(:,:) = 0.d0
       ncells(:,:) = 0.d0
       
       do i=1,nxl_full
          do j=1,nyl_full
            if (tempmask_fulldomain(i,j) .gt. 0.) then
               !accumulate values to be averaged
               global(ups%gboxx(i,j),ups%gboxy(i,j)) = global(ups%gboxx(i,j),ups%gboxy(i,j)) &
                                  + local_fulldomain(i,j)*tempmask_fulldomain(i,j)
             
               !accumulate counter that determines how many cells are being used in the average.
               !This accumulator only counts points that are included in the mask, and as such
               !avoids counting up points that are outside the 'area of interest'.        
               ncells(ups%gboxx(i,j),ups%gboxy(i,j)) = ncells(ups%gboxx(i,j),ups%gboxy(i,j)) + 1.
             end if
          enddo
       enddo
       
       !Calculate average value.
       where (ncells /= 0)            
          global = global / ncells
       elsewhere
          global(:,:) = 0.d0
       endwhere     

    end if  ! main_task

    if (allocated(local_fulldomain)) deallocate(local_fulldomain)
    if (allocated(tempmask_fulldomain)) deallocate(tempmask_fulldomain)
    if (allocated(ncells)) deallocate(ncells)
    
  end subroutine local_to_global_avg

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine local_to_global_sum(ups,local,global,mask)

    !*FD Upscale to global domain by summing local field.
    !*FD The result is an accumulated sum, not an average.
    !*FD
    !*FD Note that:
    !*FD \begin{itemize}
    !*FD \item \texttt{global} output is only valid on the main task
    !*FD \item \texttt{ups} input only needs to be valid on the main task
    !*FD \item \texttt{gboxx} and \texttt{gboxy} are the same size as \texttt{local_fulldomain}
    !*FD \item \texttt{gboxn} is the same size as \texttt{global}
    !*FD \end{itemize}

    use parallel, only : main_task, distributed_gather_var
    use nan_mod , only : NaN

    ! Arguments

    type(upscale),          intent(in)  :: ups    !*FD Upscaling indexing data.
    real(dp),dimension(:,:),intent(in)  :: local  !*FD Data on projected grid (input).
    real(dp),dimension(:,:),intent(out) :: global !*FD Data on global grid (output).
    integer,dimension(:,:),intent(in),optional :: mask !*FD Mask for upscaling

    ! Internal variables

    integer :: nxl_full,nyl_full,i,j,nxg,nyg
    real(dp),dimension(size(local,1),size(local,2)) :: tempmask

    ! values of 'local' and 'tempmask' spanning full domain (all tasks)
    real(dp),dimension(:,:), allocatable            :: local_fulldomain
    real(dp),dimension(:,:), allocatable            :: tempmask_fulldomain

    ! Beginning of code

    global = NaN

    if (present(mask)) then
       tempmask = mask
    else
       tempmask = 1.d0
    endif

    ! Gather 'local' and 'tempmask' onto main task, which is the only one that does the regridding

    call distributed_gather_var(local, local_fulldomain)
    call distributed_gather_var(tempmask, tempmask_fulldomain)

    ! Main task does regridding
    if (main_task) then

       nxl_full = size(local_fulldomain,1)
       nyl_full = size(local_fulldomain,2)
       global = 0.d0

       do i=1,nxl_full
          do j=1,nyl_full
             global(ups%gboxx(i,j),ups%gboxy(i,j)) = global(ups%gboxx(i,j),ups%gboxy(i,j)) &
                                                   + local_fulldomain(i,j)*tempmask_fulldomain(i,j)
          enddo
       enddo

    end if  ! main_task

    if (allocated(local_fulldomain)) deallocate(local_fulldomain)
    if (allocated(tempmask_fulldomain)) deallocate(tempmask_fulldomain)

  end subroutine local_to_global_sum
  
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine local_to_global_min(ups,local,global,mask)

    !*FD Upscale to global domain by finding the minimum of the local field.
    !*FD The result is an accumulated sum, not an average.
    !*FD
    !*FD Note that:
    !*FD \begin{itemize}
    !*FD \item \texttt{global} output is only valid on the main task
    !*FD \item \texttt{ups} input only needs to be valid on the main task
    !*FD \item \texttt{gboxx} and \texttt{gboxy} are the same size as \texttt{local_fulldomain}
    !*FD \item \texttt{gboxn} is the same size as \texttt{global}
    !*FD \end{itemize}

    use parallel, only : main_task, distributed_gather_var
    use nan_mod , only : NaN

    ! Arguments

    type(upscale),          intent(in)  :: ups    !*FD Upscaling indexing data.
    real(dp),dimension(:,:),intent(in)  :: local  !*FD Data on projected grid (input).
    real(dp),dimension(:,:),intent(out) :: global !*FD Data on global grid (output).
    integer,dimension(:,:),intent(in),optional :: mask !*FD Mask for upscaling

    ! Internal variables

    integer :: nxl_full,nyl_full,i,j,nxg,nyg
    real(dp),dimension(size(local,1),size(local,2)) :: tempmask

    ! values of 'local' and 'tempmask' spanning full domain (all tasks)
    real(dp),dimension(:,:), allocatable            :: local_fulldomain
    real(dp),dimension(:,:), allocatable            :: tempmask_fulldomain

    ! Beginning of code

    global = NaN

    if (present(mask)) then
       tempmask = mask
    else
       tempmask = 1
    endif

    ! Gather 'local' and 'tempmask' onto main task, which is the only one that does the regridding

    call distributed_gather_var(local, local_fulldomain)
    call distributed_gather_var(tempmask, tempmask_fulldomain)

    ! Main task does regridding
    if (main_task) then

       nxl_full = size(local_fulldomain,1)
       nyl_full = size(local_fulldomain,2)
       global(:,:) = 0.d0

       !set topography value in global cells for which the mask exists, to a very high value.
       !this should then be reduced on the next swing through the loop.
       do i=1,nxl_full
          do j=1,nyl_full    
             if (tempmask_fulldomain(i,j) > 0.) then
               global(ups%gboxx(i,j),ups%gboxy(i,j)) = huge(1.d0)
             endif
          enddo
       enddo         

       do i=1,nxl_full
          do j=1,nyl_full
             if (tempmask_fulldomain(i,j) > 0.) then
               global(ups%gboxx(i,j),ups%gboxy(i,j)) = min ( &
                 global(ups%gboxx(i,j),ups%gboxy(i,j)), &
                 local_fulldomain(i,j))
             end if
          enddo
       enddo

    end if  ! main_task

    if (allocated(local_fulldomain)) deallocate(local_fulldomain)
    if (allocated(tempmask_fulldomain)) deallocate(tempmask_fulldomain)

  end subroutine local_to_global_min  

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real(dp) function bilinear_interp(xp,yp,f)

    !*FD Performs bilinear interpolation 
    !*FD in a rectangular domain. Note that the bilinear interpolation formula is:
    !*FD  \[f_{\mathtt{x},\mathtt{y}}=(1-X')(1-Y')f_{1}+X'(1-Y')f_{2}+X'Y'f_{3}+(1-X')Y'f_{4}\]
    !*FD where $X'$ and $Y'$ are the fractional displacements of the target point within the domain.
    !*RV The value of \texttt{f} at \texttt{x,y}

    ! Argument declarations

    real(dp),             intent(in) :: xp    !*FD The fractional $x$-displacement of the target.
    real(dp),             intent(in) :: yp    !*FD The fractional $y$-displacement of the target.
    real(dp),dimension(4),intent(in) :: f     !*FD The interpolation domain;
                                              !*FD i.e. the four points surrounding the
                                              !*FD target, presented anticlockwise from bottom-
                                              !*FD left 
    ! Apply bilinear interpolation formula

    bilinear_interp=(1-xp)*(1-yp)*f(1)+xp*(1-yp)*f(2)+xp*yp*f(3)+(1-xp)*yp*f(4)

  end function bilinear_interp

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine find_ll_index(il,jl,lon,lat,lons,lats)

    !*FD Find the global gridpoint at the first corner of the box surrounding
    !*FD a given location in lat-lon space.

    use glimmer_utils

    ! Arguments

    real(dp),             intent(in)  :: lon    !*FD Longitude of location to be indexed (input)
    real(dp),             intent(in)  :: lat    !*FD Latitude of location to be indexed (input)
    real(dp),dimension(:),intent(in)  :: lats   !*FD Latitudes of global grid points 
    real(dp),dimension(:),intent(in)  :: lons   !*FD Longitudes of global grid points 
    integer,              intent(out) :: il     !*FD $x$-gridpoint index (output)
    integer,              intent(out) :: jl     !*FD $y$-gridpoint index (output)

    ! Internal variables

    integer :: nx,ny,i

    nx=size(lons) ; ny=size(lats)

    il=nx

    do i=1,nx-1
       if (lon_between(lons(i),lons(i+1),lon)) then
          il=i
          exit
       endif
    enddo

    if ((lat<lats(ny)).and.(lat>-90.0)) then
       jl=ny
       return
    endif

    if ((lat>lats(1)).and.(lat<90.0)) then
       jl=1
       return
    endif

    jl=1
    do 
       if (lat>lats(jl)) exit
       jl=jl+1
    enddo

  end subroutine find_ll_index

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine index_local_boxes (xloc, yloc, xfrac, yfrac, ggrid, proj, lgrid, lmask)

    !*FD Indexes the corners of the
    !*FD global grid box in which each local grid box sits.

    use glimmer_utils
    use glint_global_grid
    use glimmer_coordinates
    use glimmer_map_trans
    use parallel, only: main_task

    ! Arguments

    integer, dimension(:,:,:),intent(out) :: xloc,yloc   !*FD Array of indicies (see \texttt{downscale} type)
    real(dp),dimension(:,:),  intent(out) :: xfrac,yfrac !*FD Fractional off-sets of grid points
    type(global_grid),        intent(in)  :: ggrid       !*FD Global grid to be used
    type(glimmap_proj),       intent(in)  :: proj        !*FD Projection to be used
    type(coordsystem_type),   intent(in)  :: lgrid       !*FD Local grid
    integer, dimension(:,:),  intent(out) :: lmask  !*FD Mask of local cells for which interpolation is valid

    ! Internal variables

    integer :: i,j,il,jl,temp
    real(dp) :: ilon,jlat,xa,ya,xb,yb,xc,yc,xd,yd
    integer :: nx, ny, nxg, nyg, n

    if (GLC_DEBUG .and. main_task) then
       nx = lgrid%size%pt(1)
       ny = lgrid%size%pt(2)
       nxg = size(ggrid%mask,1)
       nyg = size(ggrid%mask,2)

       write(stdout,*) ' '
       write(stdout,*) 'nx,  ny =', nx, ny
       write(stdout,*) 'nxg, nyg =', nxg, nyg
       write(stdout,*) 'Indexing local boxes'
    end if

    do i=1,lgrid%size%pt(1)
       do j=1,lgrid%size%pt(2)

          ! Find out where point i,j is in lat-lon space

          call glimmap_xy_to_ll(ilon,jlat,real(i,kind=dp),real(j,kind=dp),proj,lgrid)

          ! Index that location onto the global grid

          call find_ll_index(il,jl,ilon,jlat,ggrid%lons,ggrid%lats)

          xloc(i,j,1)=il  ! This is the starting point - we now need to find
          yloc(i,j,1)=jl  ! three other points that enclose the interpolation target

          if (jlat>ggrid%lats(ggrid%ny)) then

             ! For all points except on the bottom row

             xloc(i,j,2)=il+1
             yloc(i,j,2)=jl

             xloc(i,j,3)=il+1
             yloc(i,j,3)=jl-1

             xloc(i,j,4)=il
             yloc(i,j,4)=jl-1

             call fix_bcs2d(xloc(i,j,2) ,yloc(i,j,2),ggrid%nx,ggrid%ny)
             call fix_bcs2d(xloc(i,j,3) ,yloc(i,j,3),ggrid%nx,ggrid%ny)
             call fix_bcs2d(xloc(i,j,4) ,yloc(i,j,4),ggrid%nx,ggrid%ny)

             if (jl==1) then
                temp=xloc(i,j,3)
                xloc(i,j,3)=xloc(i,j,4)
                xloc(i,j,4)=temp
             endif

          else

             ! The bottom row

             xloc(i,j,2)=il-1
             yloc(i,j,2)=jl

             xloc(i,j,3)=il-1
             yloc(i,j,3)=jl+1

             xloc(i,j,4)=il
             yloc(i,j,4)=jl+1

             call fix_bcs2d(xloc(i,j,2) ,yloc(i,j,2),ggrid%nx,ggrid%ny)
             call fix_bcs2d(xloc(i,j,3) ,yloc(i,j,3),ggrid%nx,ggrid%ny)
             call fix_bcs2d(xloc(i,j,4) ,yloc(i,j,4),ggrid%nx,ggrid%ny)

             temp=xloc(i,j,3)
             xloc(i,j,3)=xloc(i,j,4)
             xloc(i,j,4)=temp

          endif

          ! Now, find out where each of those points is on the projected
          ! grid, and calculate fractional displacements accordingly

          call glimmap_ll_to_xy(ggrid%lons(xloc(i,j,1)),ggrid%lats(yloc(i,j,1)),xa,ya,proj,lgrid)
          call glimmap_ll_to_xy(ggrid%lons(xloc(i,j,2)),ggrid%lats(yloc(i,j,2)),xb,yb,proj,lgrid)
          call glimmap_ll_to_xy(ggrid%lons(xloc(i,j,3)),ggrid%lats(yloc(i,j,3)),xc,yc,proj,lgrid)
          call glimmap_ll_to_xy(ggrid%lons(xloc(i,j,4)),ggrid%lats(yloc(i,j,4)),xd,yd,proj,lgrid)

          call calc_fractional(xfrac(i,j),yfrac(i,j),real(i,kind=dp),real(j,kind=dp), &
               xa,ya,xb,yb,xc,yc,xd,yd)

          ! If all four global gridcells surrounding this local gridcell
          ! are masked out, then mask out the local gridcell

          if ( (ggrid%mask(xloc(i,j,1),yloc(i,j,1)) == 0) .and.  &
               (ggrid%mask(xloc(i,j,2),yloc(i,j,2)) == 0) .and.  &
               (ggrid%mask(xloc(i,j,3),yloc(i,j,3)) == 0) .and.  &
               (ggrid%mask(xloc(i,j,4),yloc(i,j,4)) == 0) ) then
             lmask(i,j) = 0
          else
             lmask(i,j) = 1
          endif

       enddo
    enddo

    ! This output should work for any number of tasks, since the variables here span the
    ! full local domain. Because of where it falls in the call chain, only the main task
    ! should reach this point, but we check that anyway to avoid problems.

    !WHL - Commenting out these diagnostics for now
    !TODO - Modify this code for multiple tasks

!!    if (GLC_DEBUG .and. main_task) then

!!       write(stdout,*) ' '
!!       write(stdout,*) 'Mask in neighborhood of i, j = ', itest_local, jtest_local
!!       do j = jtest_local-1, jtest_local+1
!!          write(stdout,*) lmask(itest_local-1:itest_local+1,j)
!!       enddo
     
!!       write(stdout,*) ' '
!!       write(stdout,*) 'Global mask near Greenland'
!!       do j = 1, 20
!!          ! This is hardwired for a particular global grid 
!!          write(stdout,150) ggrid%mask(nxg-29:nxg,j)   
!!       enddo

!!       write(stdout,*) ' '
!!       write(stdout,*) 'Local mask'
!!       do j = ny, 1, -1
!!          write(stdout,200) lmask(1:nx,j)
!!       enddo

!!  150  format(30i2)
!!  200  format(76i2)
!!    end if

  end subroutine index_local_boxes

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine calc_grid_angle(downs,proj,lgrid)

    !*FD Calculates the angle the projected 
    !*FD grid makes with north at each point and stores the cos 
    !*FD and sin of that angle in the relevant arrays in \texttt{proj}.

    use glimmer_coordinates
    use glimmer_map_trans

    type(downscale),intent(inout) :: downs !*FD The projection to be used
    type(glimmap_proj),intent(in) :: proj
    type(coordsystem_type),intent(in) :: lgrid

    integer :: i,j
    real(dp) :: latn,lonn,lats,lons,lat,lon,dlat,dlon,temp

    do i=1,lgrid%size%pt(1)

       ! Main, central block

       do j=2,lgrid%size%pt(2)-1
          call glimmap_xy_to_ll(lonn,latn,real(i,kind=dp),real(j+1,kind=dp),proj,lgrid)
          call glimmap_xy_to_ll(lon,lat,real(i,kind=dp),real(j,kind=dp),proj,lgrid)
          call glimmap_xy_to_ll(lons,lats,real(i,kind=dp),real(j-1,kind=dp),proj,lgrid)
          dlat=latn-lats
          dlon=lonn-lons
          if (dlon<-90) dlon=dlon+360
          temp=atan(dlon/dlat)
          downs%sintheta(i,j)=sin(temp)
          downs%costheta(i,j)=cos(temp)
       enddo

       ! bottom row

       call glimmap_xy_to_ll(lonn,latn,real(i,kind=dp),real(2,kind=dp),proj,lgrid)
       call glimmap_xy_to_ll(lon,lat,real(i,kind=dp),real(1,kind=dp),proj,lgrid)
       dlat=latn-lat
       dlon=lonn-lon
       if (dlon<-90) dlon=dlon+360
       temp=atan(dlon/dlat)
       downs%sintheta(i,1)=sin(temp)
       downs%costheta(i,1)=cos(temp)

       ! top row

       call glimmap_xy_to_ll(lon,lat,real(i,kind=dp),real(lgrid%size%pt(2),kind=dp),proj,lgrid)
       call glimmap_xy_to_ll(lons,lats,real(i,kind=dp),real(lgrid%size%pt(2)-1,kind=dp),proj,lgrid)
       dlat=lat-lats
       dlon=lon-lons
       if (dlon<-90) dlon=dlon+360
       temp=atan(dlon/dlat)
       downs%sintheta(i,lgrid%size%pt(2))=sin(temp)
       downs%costheta(i,lgrid%size%pt(2))=cos(temp)

    enddo

  end subroutine calc_grid_angle

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine new_upscale(ups,ggrid,proj,mask,lgrid)

    use glint_global_grid
    use glimmer_log
    use glimmer_map_trans
    use glimmer_coordinates

    !*FD Compiles an index of which global grid box contains a given
    !*FD grid box on the projected grid, and sets derived type \texttt{ups}
    !*FD accordingly.

    ! Arguments

    type(upscale),         intent(out) :: ups        !*FD Upscaling type to be set
    type(global_grid),     intent(in)  :: ggrid      !*FD Global grid to be used
    type(glimmap_proj),    intent(in)  :: proj       !*FD Projection being used
    integer,dimension(:,:),intent(in)  :: mask       !*FD Upscaling mask to be used
    type(coordsystem_type),intent(in)  :: lgrid      !*FD local grid

    ! Internal variables

    integer  :: i,j,ii,jj,nx,ny,gnx,gny
    real(dp) :: plon,plat

    ! Beginning of code

    if (associated(ups%gboxx)) deallocate(ups%gboxx)
    if (associated(ups%gboxy)) deallocate(ups%gboxy)
    if (associated(ups%gboxn)) deallocate(ups%gboxn)

    allocate(ups%gboxx(lgrid%size%pt(1),lgrid%size%pt(2)))
    allocate(ups%gboxy(lgrid%size%pt(1),lgrid%size%pt(2)))     
    allocate(ups%gboxn(ggrid%nx,ggrid%ny))

    gnx=ggrid%nx ; gny=ggrid%ny
    nx =lgrid%size%pt(1) ; ny =lgrid%size%pt(2)

    ups%gboxx=0 ; ups%gboxy=0

    do i=1,nx
       do j=1,ny
          call glimmap_xy_to_ll(plon,plat,real(i,kind=dp),real(j,kind=dp),proj,lgrid)
          ii=1 ; jj=1
          do
             ups%gboxx(i,j)=ii
             if (ii>gnx) then
                call write_log('global index failure',GM_FATAL,__FILE__,__LINE__)
             endif
             if (lon_between(ggrid%lon_bound(ii),ggrid%lon_bound(ii+1),plon)) exit
             ii=ii+1
          enddo

          jj=1

          do
             ups%gboxy(i,j)=jj
             if (jj>gny) then
                call write_log('global index failure',GM_FATAL,__FILE__,__LINE__)
             endif
             if ((ggrid%lat_bound(jj)>=plat).and.(plat>ggrid%lat_bound(jj+1))) exit
             jj=jj+1
          enddo

       enddo
    enddo

    ups%gboxn=0

    do i=1,nx
       do j=1,ny
          ups%gboxn(ups%gboxx(i,j),ups%gboxy(i,j))=ups%gboxn(ups%gboxx(i,j),ups%gboxy(i,j))+mask(i,j)
       enddo
    enddo

    ups%set=.true.

  end subroutine new_upscale

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine copy_upscale(in,out)

    use glimmer_log

    type(upscale),intent(in)  :: in
    type(upscale),intent(out) :: out

    if (.not.in%set) then
       call write_log('Attempt to copy un-initialised upscale type',GM_FATAL,&
            __FILE__,__LINE__)
    endif

    if (associated(out%gboxx)) deallocate(out%gboxx)
    if (associated(out%gboxy)) deallocate(out%gboxy)
    if (associated(out%gboxn)) deallocate(out%gboxn)

    allocate(out%gboxx(size(in%gboxx,1),size(in%gboxx,2)))
    allocate(out%gboxy(size(in%gboxy,1),size(in%gboxy,2)))
    allocate(out%gboxn(size(in%gboxn,1),size(in%gboxn,2)))

    out%gboxx=in%gboxx
    out%gboxy=in%gboxy
    out%gboxn=in%gboxn

    out%set=.true.

  end subroutine copy_upscale

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  logical function lon_between(a,b,x)

    !*FD Checks to see whether a 
    !*FD longitudinal coordinate is between two bounds,
    !*FD taking into account the periodic boundary conditions.
    !*RV Returns \texttt{.true.} if $\mathtt{x}\geq \mathtt{a}$ and $\mathtt{x}<\mathtt{b}$.

    ! Arguments

    real(dp),intent(in) :: a  !*FD Lower bound on interval for checking
    real(dp),intent(in) :: b  !*FD Upper bound on interval for checking
    real(dp),intent(in) :: x  !*FD Test value (degrees)

    ! Internal variables

    real(dp) :: ta,tb

    ! Beginning of code

    if (a<b) then
       lon_between=((x>=a).and.(x<b))
    else
       if (x<a) then
          ta=a-360.0
          tb=b
       else 
          ta=a
          tb=b+360.0
       endif
       lon_between=((x>=ta).and.(x<tb))
    endif

  end function lon_between

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine calc_fractional(x,y,xp,yp,xa,ya,xb,yb,xc,yc,xd,yd)

    !*FD Performs a coordinate transformation to locate the point
    !*FD $(X',Y')$ fractionally within an arbitrary quadrilateral, 
    !*FD defined by the points $(x_A,y_A)$, $(x_B,y_B)$, 
    !*FD $(x_C,y_C)$ and $(x_D,y_D)$, which are ordered 
    !*FD anticlockwise.

    real(dp),intent(out) :: x !*FD The fractional $x$ location.
    real(dp),intent(out) :: y !*FD The fractional $y$ location.
    real(dp),intent(in)  :: xp,yp,xa,ya,xb,yb,xc,yc,xd,yd

    real(dp) :: a,b,c
    real(dp),parameter :: small=1d-8

    a=(yb-ya)*(xc-xd)-(yc-yd)*(xb-xa)

    b=xp*(yc-yd)-yp*(xc-xd) &
         +xd*(yb-ya)-yd*(xb-xa) &
         -xp*(yb-ya)+yp*(xb-xa) &
         -xa*(yc-yd)+ya*(xc-xd) 

    c=xp*(yd-ya)+yp*(xa-xd)+ya*xd-xa*yd

    if (abs(a)>small) then
       x=(-b-sqrt(b**2-4*a*c))/(2*a)
    else
       x=-c/b
    endif

    y=(yp-ya-x*(yb-ya))/(yd+x*(yc-yd-yb+ya)-ya)

    if (GLC_DEBUG) then
! Could use the following code if points are degenerate (a=b, c=d, etc.)
!       if (abs(a) > small) then
!          x=(-b-sqrt(b**2-4*a*c))/(2*a)
!       elseif (abs(b) > small) then
!          x=-c/b
!       else
!          x=0.d0
!       endif
!
!       if (abs(yd+x*(yc-yd-yb+ya)-ya) > small) then
!          y=(yp-ya-x*(yb-ya))/(yd+x*(yc-yd-yb+ya)-ya)
!       else
!          y=0.d0
!       endif
    end if

  end subroutine calc_fractional

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module glint_interp
