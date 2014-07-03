!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glide_mask.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

#include "glide_mask.inc"

module glide_mask

    ! masking ice thicknesses

    use glimmer_global, only : dp, sp
    use nan_mod, only : NaN

    implicit none

contains

!TODO - Remove iarea, ivol, and exec_serial?
!       If iarea and ivol are desired, they can be computed elsewhere.

!TODO - This subroutine is called from glissade_velo_driver with stagthck and stagtopg
!       as input arguments.  Might be safer to have a difference mask subroutine to
!       compute the staggered mask.

  subroutine glide_set_mask(numerics, thck, topg, ewn, nsn, eus, mask, iarea, ivol, exec_serial)

    use parallel
    use glide_types
    use glimmer_physcon, only : rhoi, rhoo
    implicit none

    type(glide_numerics), intent(in) :: numerics !Numerical parameters structure
    real(dp), dimension(:,:), intent(in) :: thck !Ice thickness
    real(dp), dimension(:,:), intent(in) :: topg !Bedrock topography (not lower surface!)
    integer, intent(in) :: ewn, nsn !Grid size
    real(sp), intent(in) :: eus !Sea level
    integer, dimension(:,:), intent(inout) :: mask !Output mask
    real(dp), intent(inout), optional :: ivol, iarea !Area and volume of ice
    logical, optional :: exec_serial  !JEFF If executing in serial in MPI program.

    ! local variables
    integer ew,ns
    real(dp), parameter :: con = - rhoi / rhoo
    logical :: exec_serial_flag

!TODO - This array may not be needed, at least in parallel.

    ! Create an array to "fake" the boundaries of the mask so that boundary
    ! finding can work even on the boundaries of the real mask.

    integer, dimension(0:ewn+1,0:nsn+1) :: maskWithBounds;

!TODO - What is the exec_serial option?  Is it still needed?
    !JEFF Handle exec_serial optional parameter
    if ( present(exec_serial) ) then
       exec_serial_flag = exec_serial
    else
       ! Default to off
       exec_serial_flag = .FALSE.
    endif

    mask = 0

    if (present(iarea)) iarea = 0.d0
    if (present(ivol)) ivol = 0.d0

!TODO - This mask is confusing.  Wondering if we should replace it by a series of logical masks.

! Would need the following:
! glide_mask_has_ice = 1
! glide_mask_thin_ice = 3
! glide_mask_ocean = 4 (below sea level, with or without ice)
! glide_mask_land = 8 (complement of glide_mask_ocean)
! glide_mask_grounding_line = 16 (could define in terms of margin and has ice?)
! glide_mask_margin = 32 (has_ice + at least one neighbor with no ice)
! glide_mask_dirichlet_bc = 64
! glide_mask_comp_domain_bnd = 128 (no longer needed with new global BC?)
! glide_no_ice (complement of glide_has_ice)
! glide_is_thin
! glide_is_ocean (ocean + no_ice; change to glide_ocean_icefree or remove?)
! glide_is_land (land + no_ice; change to glide_land_icefree or remove?)
! glide_is_ground (land + has_ice)
! glide_is_float (ocean + has_ice)
! glide_is_grounding_line (just inside or just outside? Used only in glide_ground)
! glide_is_margin
! glide_is_land_margin (margin + land + has_ice)
! glide_is_calving (margin + ocean + has_ice; change the name to is_marine_margin?)
! glide_is_marine_ice_edge (margin + (float or GL); may not be needed)
! glide_is_dirichlet_boundary 
! glide_is_comp_domain_bnd (may not be needed with new global BC?)
! 
!TODO - Even if we keep the present structure, could change glide_is_land to glide_icefree_land,
!                                                           glide_is_ocean to glide_icefree_ocean
!       Could get by with fewer masks in the code by removing some combinations
!       Could remove *BITS

!TODO - Combine the following into one do loop with ifs?
!       Probably should loop over locally owned cells only, then do a halo update.

    !Identify points with any ice
    where (thck > 0.d0)
        mask = ior(mask, GLIDE_MASK_HAS_ICE)  ! GLIDE_MASK_HAS_ICE = 1; see glide_mask.inc
    endwhere

    !Identify points where the ice is below the ice dynamics limit
    where (thck > 0.d0 .and. thck < numerics%thklim)
        mask = ior(mask, GLIDE_MASK_THIN_ICE)  ! GLIDE_MASK_THIN_ICE = 3
    endwhere

    !Identify points where the ice is floating or where there is open ocean
    where (topg - eus < con * thck)
        mask = ior(mask, GLIDE_MASK_OCEAN)   ! GLIDE_MASK_OCEAN = 8
    elsewhere
        mask = ior(mask, GLIDE_MASK_LAND)    ! GLIDE_MASK_LAND = 4
    endwhere

    if (present(iarea) .and. present(ivol)) then
        call get_area_vol(thck, numerics%dew, numerics%dns, numerics%thklim, iarea, ivol, exec_serial_flag)
    end if

!TODO - The following could be accomplished by a halo call for 'mask' with appropriate global BC.
!       
    maskWithBounds = 0
    maskWithBounds(1:ewn, 1:nsn) = MASK
    maskWithBounds(0,1:nsn) = mask(1,:)
    maskWithBounds(1:ewn,0) = mask(:,1)
    maskWithBounds(ewn+1,1:nsn) = mask(ewn,:)
    maskWithBounds(1:ewn,nsn+1) = mask(:,nsn)
    maskWithBounds(0,0) = mask(1,1)
    maskWithBounds(ewn+1,nsn+1) = mask(ewn,nsn)
    maskWithBounds(0,nsn+1) = mask(1,nsn)
    maskWithBounds(ewn+1,0) = mask(ewn,1)

    ! finding boundaries

!TODO - For parallel glissade code, this loop should be only over locally owned scalars.
!       If halo cells are present, maskWithBounds array may not be needed; can replace with mask array.
!TODO - Not sure what happens here when we're computing a mask on the velocity grid.

    do ns = 1,nsn
       do ew = 1,ewn
          !Find the grounding line
          if (GLIDE_IS_GROUND(MASK(ew,ns))) then    ! land + has_ice
             if (GLIDE_IS_FLOAT(maskWithBounds(ew-1,ns)) .or. &
                  GLIDE_IS_FLOAT(maskWithBounds(ew+1,ns)) .or. &
                  GLIDE_IS_FLOAT(maskWithBounds(ew,ns-1)) .or. & 
                  GLIDE_IS_FLOAT(maskWithBounds(ew,ns+1))) then
                MASK(ew,ns) = ior(MASK(ew,ns),GLIDE_MASK_GROUNDING_LINE)
             end if
          end if

          ! Ice margin
          ! *tb* A point is now masked even if it touches the ocean on one corner.
          if ( GLIDE_HAS_ICE(mask(ew, ns)) .and. &
              (GLIDE_NO_ICE(maskWithBounds(ew-1,ns))   .or. GLIDE_NO_ICE(maskWithBounds(ew+1,ns))   .or. &
               GLIDE_NO_ICE(maskWithBounds(ew,ns-1))   .or. GLIDE_NO_ICE(maskWithBounds(ew,ns+1))   .or. &
               GLIDE_NO_ICE(maskWithBounds(ew-1,ns-1)) .or. GLIDE_NO_ICE(maskWithBounds(ew-1,ns+1)) .or. &
               GLIDE_NO_ICE(maskWithBounds(ew+1,ns-1)) .or. GLIDE_NO_ICE(maskWithBounds(ew+1,ns+1)))) then
             MASK(ew,ns) = ior(MASK(ew,ns),GLIDE_MASK_MARGIN)
          end if

!TODO - Not sure if this will be needed when global boundaries are handled correctly.
!       The GLIDE_MASK_COMP_DOMAIN_BND condition is currently used in glam_strs2.F90.

         !Mark domain boundaries
         !if (ns == 1 .or. ns == nsn .or. ew == 1 .or. ew == ewn) then
         if (parallel_boundary(ew,ewn,ns,nsn)) then
! SFP: commenting out for now, while trying to get periodic bcs working
!            mask(ew, ns) = ior(mask(ew, ns), GLIDE_MASK_COMP_DOMAIN_BND)
         end if
       end do
    end do

!TODO - This halo update should be moved to a higher level.

    !JEFF Don't call halo update if running in serial mode
    if (.NOT. exec_serial_flag) then
       call parallel_halo(mask)
    endif

  end subroutine glide_set_mask

  subroutine augment_kinbc_mask(mask, kinbcmask)

    !TODO adding the kinematic bc to the unstaggered mask no longer needs to be supported.  That functionality can be removed.
    !*FD Augments the Glide mask with the location of kinematic (dirichlet) boundary
    !*FD conditions.  These locations cannot be determined by the model a priori, and
    !*FD must be specified through a field in a NetCDF file.
    integer, dimension(:,:), target :: mask
    integer, dimension(:,:) :: kinbcmask

    integer, dimension(:,:), pointer :: maskp

    !Because the kinematic boundary conditions are specified on the staggered grid,
    !there may be a size mismatch here depending on whether we are computing a mask
    !for the staggered grid.
    if (size(mask, 1) /= size(kinbcmask, 1)) then
        maskp => mask(1:size(mask,1) - 1, 1:size(mask,2) - 1)
    else
        maskp => mask
    end if

    where (kinbcmask /= 0)
        maskp = ior(maskp, GLIDE_MASK_DIRICHLET_BC)
    endwhere
  end subroutine augment_kinbc_mask

  subroutine get_area_vol(thck, dew, dns, thklim, iarea, ivol, exec_serial)
    use parallel
    implicit none
    real(dp), dimension(:,:) :: thck
    real(dp) :: dew, dns, thklim
    real(dp) :: iarea, ivol, sum(2)
    logical :: exec_serial

    integer :: i,j

    do i = 1+lhalo, size(thck,1)-uhalo
        do j = 1+lhalo, size(thck,2)-uhalo
            if (thck(i,j) > thklim ) then
                iarea = iarea + 1
                ivol = ivol + thck(i,j)
            end if
        end do
    end do

    iarea = iarea  * dew * dns
    ivol = ivol * dew * dns
    
    if (.NOT. exec_serial) then
       sum(1) = iarea
       sum(2) = ivol
       call global_sum(sum)
       iarea = sum(1)
       ivol  = sum(2)
    endif

  end subroutine get_area_vol
 
  subroutine calc_iareaf_iareag(dew, dns, iarea, mask, iareaf, iareag, exec_serial)
    
    use parallel
    !TODO - remove iarea from the call since it is not used

    implicit none
    real(dp), intent(in) :: dew, dns
    real(dp), intent(in) :: iarea
    real(dp), intent(out) :: iareaf, iareag
    integer, dimension(:,:), intent(in) :: mask 
    logical, optional :: exec_serial  ! If executing in serial in MPI program.

    integer :: i,j
    logical :: exec_serial_flag
    real(dp) :: sum(2)
 
    !TODO - Is this exec_serial option needed?
    ! Handle exec_serial optional parameter
    if ( present(exec_serial) ) then
      exec_serial_flag = exec_serial
    else
      ! Default to off
      exec_serial_flag = .FALSE.
    endif

    iareaf = 0.d0
    iareag = 0.d0 

    !loop over locally owned scalars
    do j = 1+lhalo, size(mask,2)-uhalo
      do i = 1+lhalo, size(mask,1)-uhalo
        if (GLIDE_IS_FLOAT(mask(i,j))) then
          iareaf = iareaf + dew * dns
        else if(GLIDE_IS_GROUND_OR_GNDLINE(mask(i,j))) then
          iareag = iareag + dew * dns
        end if
      end do
    end do

    if (.NOT. exec_serial_flag) then
       sum(1) = iareaf
       sum(2) = iareag
       call global_sum(sum)
       iareaf = sum(1)
       iareag = sum(2)
    endif

  end subroutine calc_iareaf_iareag

    subroutine glide_marine_margin_normal(thck, mask, marine_bc_normal, exec_serial)

      !TODO - Steve thinks this is a PBJ routine that could be removed.

      use parallel
        use glimmer_physcon, only:pi
!!        use glimmer_horiz_bcs, only: horiz_bcs_unstag_scalar
        implicit none
        !*FD This subroutine derives from the given mask the normal to an ice shelf
        !*FD each point on the marine margin.
        real(dp), dimension(:,:), intent(in) :: thck
        integer, dimension(:,:), intent(in) :: mask
        real(dp), dimension(:,:), intent(out) :: marine_bc_normal
        logical, optional :: exec_serial  !JEFF If executing in serial in MPI program.

        integer :: i, j, dx, dy, k
        logical :: exec_serial_flag

        real(dp), dimension(size(thck,1), size(thck,2)) :: direction_x, direction_y
        
        real(dp), dimension(-1:1, -1:1) :: angle_lookup

	    !JEFF Handle exec_serial optional parameter
	    if ( present(exec_serial) ) then
	       exec_serial_flag = exec_serial
	    else
	       ! Default to off
	       exec_serial_flag = .FALSE.
	    endif

                !direction_y =    -1       0       1        !direction_x = 
        angle_lookup(-1, :) = (/ 3*pi/4,   pi/2,   pi/4 /)  !-1
        angle_lookup( 0, :) = (/   pi,     0D0,  2*pi   /)  ! 0
        angle_lookup( 1, :) = (/ 5*pi/4, 3*pi/2, 7*pi/4 /)  ! 1
        call upwind_from_mask(mask, direction_x, direction_y, exec_serial_flag)

        !Set up a thickness variable with "ghost cells" so that we don't go out
        !of bounds with the vectorized operation below
        !thckWithBounds(1:size(thck,1), 1:size(thck,2)) = thck
        !thckWithBounds(:,0) = thckWithBounds(:,1)
        !thckWithBounds(0,:) = thckWithBounds(1,:)
        !thckWithBounds(size(thck,1)+1,:) = thckWithBounds(size(thck,1),:)
        !thckWithBounds(:,size(thck,2)+1) = thckWithBounds(:,size(thck,2))
        do i = 1, size(mask, 1)
            do j = 1, size(mask, 2)
                if (GLIDE_IS_CALVING(mask(i,j))) then
                    dx = int(direction_x(i,j))
                    dy = int(direction_y(i,j))
                    if (dx == 0 .and. dy == 0) then
                        write(*,*)"A shelf front point has been identified at:"
                        write(*,*)"x = ",i
                        write(*,*)"y = ",j
                        write(*,*)"But neither x nor y derivatives have been marked as upwinded."
                        write(*,*)"This should never happen, if this error appears it is a bug"
                        write(*,*)"and should be reported."
                        write(*,*)"The mask around this point follows:"
                        write(*,*)"--------------------------"

                        !Write a header row with a * in the column corresponding to the center
                        do k = -4, 4
                            if (k==0) then
                                write(*,"(A)",advance="no")"           *"
                            else if (i+k > 0 .and. i+k <= size(mask,1)) then
                                write(*,"(A)",advance="no")"            "
                            end if
                        end do 
                        write(*,*)

                        do k=4, -4, -1
                            if (j+k > 0 .and. j+k <= size(mask, 2)) then
                                if (k == 0) then
                                    write(*,*) "*", mask(max(1,i-4):min(size(mask,1),i+4),j+k)
                                else
                                    write(*,*) " ", mask(max(1,i-4):min(size(mask,1),i+4),j+k)
                                end  if
                            end if
                        end do
                        write(*,*)"--------------------------"
                        write(*,*)"Have a nice day!"
                        !stop
                    end if
                    marine_bc_normal(i,j) = angle_lookup(dx, dy) 
                    !marine_bc_normal(i,j) = calc_normal_45deg(thckWithBounds(i-1:i+1,j-1:j+1))
                else
                    marine_bc_normal(i,j) = NaN
                end if
            end do
        end do
        if (.NOT. exec_serial_flag) then
           call parallel_halo(marine_bc_normal)
!           call horiz_bcs_unstag_scalar(marine_bc_normal)
        endif
    end subroutine

    !TODO - This is a PBJ function and could be removed.
    function calc_normal_45deg(thck3x3)
        use glimmer_physcon, only: pi
        
        !*FD Computes the angle of the normal vector, in radians, for the given
        !*FD 3x3 segment of ice geometry.
        !*FD The normal is given in increments of 45 degrees (no nicer
        !*FD interpolation is currently done)
        !*FD This is based on the Payne and Price GLAM code, if/when this is
        !*FD integrated into CISM it should probably be refactored to use this.
        real(dp), dimension(3,3) :: thck3x3

        real(dp) :: calc_normal_45deg
         
        real(dp), dimension(3,3) :: mask, maskcorners
        real(dp), dimension(3,3) :: thckmask
        real(dp), dimension(3) :: testvect
        real(dp) :: phi, deg2rad
        integer :: loc_latbc

        deg2rad = pi / 180.0d0
        loc_latbc = 0
        phi = 0.d0
        mask(:,1) = (/ 0.0d0, 180.0d0, 0.0d0 /)
        mask(:,2) = (/ 270.0d0, 0.0d0, 90.0d0 /)
        mask(:,3) = (/ 0.0d0, 360.0d0, 0.0d0 /)
        maskcorners(:,1) = (/ 225.0d0, 0.0d0, 135.0d0 /)
        maskcorners(:,2) = (/ 0.0d0, 0.0d0, 0.0d0 /)
        maskcorners(:,3) = (/ 315.0d0, 0.0d0, 45.0d0 /)

        ! specify new value of 'loc' vector such that fwd/bwd diffs. are set up correctly in sparse matrix
        ! when function 'fillsprsebndy' is called. Also, specify appropriate values for the vectors 'normal'
        ! and 'fwdorbwd', which specify the orientation of the boundary normal and the direction of forward or
        ! backward differencing to be done in the lateral boundary condition functions 'normhorizmainbc_lat'
        ! and 'crosshorizmainbc_lat'

        ! following is algorithm for calculating boundary normal at 45 deg. increments, based on arbitray
        ! boundary shape

        where( thck3x3  /=  0.0d0 )
            thckmask = 0.0_dp
        elsewhere( thck3x3 == 0.0d0 )
            thckmask = 1.0d0
        endwhere

        testvect = sum( thckmask * mask, 1 )

        !if( up == 3 )then ! temporary code for debugging
        !  do i = 3,1,-1
        !  print *, 'thck = ', thck(:,i)
        !  end do
        !  print *, ' '
        !
        !  do i = 3,1,-1
        !      print *, 'thckmask = ', thckmask(:,i)
        !  end do
        !  print *, ' '
        !
        !  print *, 'testvect =  ', testvect
        !  print *, ' '
        !end if

        ! calculate the angle of the normal in cart. (x,y) system w/ 0 deg. at 12 O'clock, 90 deg. at 3 O'clock, etc.
        if( sum( sum( thckmask, 1 ) ) == 1.0d0 )then
            phi = sum( sum( thckmask * maskcorners, 1 ) )
        else
            if( any( testvect == 360.0d0 ) )then
                if( sum( testvect ) == 450.0d0 )then
                    phi = 45.0d0
                elseif( sum( testvect ) == 630.0d0 )then
                    phi = 315.0d0
                else
                    phi = 0.0d0
                end if
            elseif( all( testvect  /=  360 ) )then
                phi = sum( testvect ) / sum( testvect/testvect, testvect  /=  0.0d0 )
            end if
        end if

        calc_normal_45deg = deg2rad * phi
        
        !Tim's Note: This appears to actually compute 0 at 6 O'clock according
        !to Glimmer's coordinate system.  90 deg. is still 3 O'clock.
        !I'm going to correct for this here rather than dig through the code
        !above
        !(TODO: correct it in the code above!)
        calc_normal_45deg = pi - calc_normal_45deg 
        if (calc_normal_45deg < 0) calc_normal_45deg = calc_normal_45deg + 2*pi

    end function

!TODO - This subroutine may not be needed.
!       It is not currently called from anywhere.

    !Fills a field of differencing directions suitable to give a field
    !derivative routine.  Uses centered differencing everywhere except for the
    !marine ice margin, where upwinding and downwinding is used to avoid
    !differencing across the boundary.

    subroutine upwind_from_mask(mask, direction_x, direction_y, exec_serial)
      use parallel
!!        use glimmer_horiz_bcs, only: horiz_bcs_unstag_scalar
        integer, dimension(:,:), intent(in) :: mask
        double precision, dimension(:,:), intent(out) :: direction_x, direction_y
        logical, optional :: exec_serial  !JEFF If executing in serial in MPI program.

        integer :: i,j
        logical :: exec_serial_flag

	    !JEFF Handle exec_serial optional parameter
	    if ( present(exec_serial) ) then
	       exec_serial_flag = exec_serial
	    else
	       ! Default to off
	       exec_serial_flag = .FALSE.
	    endif

        direction_x = 0
        direction_y = 0

        !Detect locations of the marine margin
        do i = 1, size(mask,1)
            do j = 1, size(mask,2)
                if (GLIDE_IS_CALVING(mask(i,j))) then
                    !Detect whether we need to upwind or downwind in the Y
                    !direction
                    if (i > 1) then
                        if (.not. GLIDE_HAS_ICE(mask(i-1,j))) then
                            direction_x(i,j) = 1
                        end if
                    end if

                    if (i < size(mask, 1)) then
                        if (.not. GLIDE_HAS_ICE(mask(i+1,j))) then
                            direction_x(i,j) = -1
                        end if
                    end if

                    !Detect whether we need to upwind or downwind in the X
                    !direction
                    if (j > 1) then
                        if (.not. GLIDE_HAS_ICE(mask(i,j-1))) then
                            direction_y(i,j) = 1
                        end if
                    end if
                    
                    if (j < size(mask, 2)) then
                        if (.not. GLIDE_HAS_ICE(mask(i,j+1))) then
                            direction_y(i,j) = -1
                        end if
                    end if

                    !If we are at a point that is "interior" to two other boundary points, 
                    !such as the lower right of:
                    !o b i
                    !b b i
                    !(o = ocean, b = boundary, i = interior), then we will not detect the need
                    !to upwind or downwind.  However, we still should for consistency with other
                    !mask points (in some cases, not doing so can lead to a singular calculation
                    !at the marine ice front)
                    !
                    !We can think of this operation as avoiding calving points where there is 
                    !a non-calving point to upwind into.
                    !
                    !TODO: We need a better way to detect interior points.  Right now I am just using
                    !points that are floating, and that works, but this doesn't work for two reasons:
                    !1. Boundary points are also floating
                    !2. Could fail for a very thin ice shelf
                    if (int(direction_x(i,j)) == 0 .and. int(direction_y(i,j)) == 0 .and. &
                        i > 1 .and. j > 1 .and. i < size(mask, 1) .and. j < size(mask, 2)) then
                        if (.not. GLIDE_HAS_ICE(mask(i-1, j-1))) then
                            direction_x(i,j) = 1
                            direction_y(i,j) = 1
                        else if (.not. GLIDE_HAS_ICE(mask(i-1, j+1))) then
                            direction_x(i,j) = 1
                            direction_y(i,j) = -1
                        else if (.not. GLIDE_HAS_ICE(mask(i+1, j-1))) then
                            direction_x(i,j) = -1
                            direction_y(i,j) = 1
                        else if (.not. GLIDE_HAS_ICE(mask(i+1, j+1))) then
                            direction_x(i,j) = -1
                            direction_y(i,j) = -1
                        end if
                    end if
                end if
            end do
        end do

        if (.NOT. exec_serial_flag) then
            call parallel_halo(direction_x)
!            call horiz_bcs_unstag_scalar(direction_x)
            call parallel_halo(direction_y)
!            call horiz_bcs_unstag_scalar(direction_y)
        endif
    end subroutine upwind_from_mask

end module glide_mask
