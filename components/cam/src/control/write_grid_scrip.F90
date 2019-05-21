module write_grid_scrip_mod
   use netcdf
   use shr_kind_mod,    only: r8 => shr_kind_r8
   use cam_abortutils,  only: endrun
   implicit none
   
   public :: write_grid_scrip

   private :: check

contains
!===================================================================================================
!===================================================================================================
subroutine check(status)
   use cam_logfile,    only: iulog
   integer, intent ( in) :: status
   if(status /= NF90_NOERR) then 
      write(iulog,*), trim(nf90_strerror(status))
      write(iulog,*) "ERROR - write_grid_scrip - status: ",status 
      call endrun("ERROR - write_grid_scrip" )
   end if
end subroutine check 
!===================================================================================================
!===================================================================================================
subroutine write_grid_scrip()
   use spmd_utils,         only: masterproc
   use dimensions_mod,     only: nelem, nelemd, nelemdmax, ne, np, npsq
   use dyn_grid,           only: compute_global_area, compute_global_coords, fv_nphys
   use phys_grid,          only: ngcols
   !----------------------------------------------------------------------------
   character(len=100) :: filename
   integer, parameter :: nvar = 7
   integer            :: ncid
   integer            :: varid(nvar)
   integer            :: dimid_column
   integer            :: dimid_corner
   integer            :: dimid_rank
   integer            :: dimid_2d(2)
   integer            :: icol, c
   integer            :: lon_diff
   real(r8), pointer             :: area(:)
   integer,  dimension(ngcols)   :: mask
   real(r8), dimension(ngcols)   :: center_lat_rad 
   real(r8), dimension(ngcols)   :: center_lon_rad 
   real(r8), dimension(ngcols)   :: center_lat_deg 
   real(r8), dimension(ngcols)   :: center_lon_deg 
   real(r8), dimension(ngcols,5) :: corner_lat_deg 
   real(r8), dimension(ngcols,5) :: corner_lon_deg 
   real(r8), dimension(5,ngcols) :: corner_lat_deg_alt 
   real(r8), dimension(5,ngcols) :: corner_lon_deg_alt 
   !----------------------------------------------------------------------------

   mask(:) = 1

   ! Retreive grid properties
   allocate(area(ngcols))
   call compute_global_area(area)

   call compute_global_coords(center_lat_rad, center_lon_rad, &
                              center_lat_deg, center_lon_deg, &
                              corner_lat_deg, corner_lon_deg )
   
   ! Add redundant corner
   corner_lat_deg(:,5) = corner_lat_deg(:,1)
   corner_lon_deg(:,5) = corner_lon_deg(:,1)

   ! reshape for proper dimension ordering
   do icol = 1,ngcols
      corner_lat_deg_alt(:,icol) = corner_lat_deg(icol,:)
      corner_lon_deg_alt(:,icol) = corner_lon_deg(icol,:)
      do c = 1,5
         lon_diff = center_lon_deg(icol) - corner_lon_deg_alt(c,icol)
         if (lon_diff > 180.) then
            corner_lon_deg_alt(c,icol) = corner_lon_deg_alt(c,icol) + 360.
         end if
         if (lon_diff < -180.) then
            corner_lon_deg_alt(c,icol) = corner_lon_deg_alt(c,icol) - 360.
         end if
      end do 
   end do
   !----------------------------------------------------------------------------
   ! Only use masterproc to write the file
   if (masterproc) then

      write(filename,"(A2,I0,A2,I0,A9)")"ne",ne,"pg",fv_nphys,"_scrip.nc"

      ! Create the file
      call check( nf90_create(filename, NF90_CLOBBER, ncid) )

      ! Define the dimensions
      call check( nf90_def_dim(ncid, "grid_size",    ngcols, dimid_column) )
      call check( nf90_def_dim(ncid, "grid_corners", 5,      dimid_corner) )
      call check( nf90_def_dim(ncid, "grid_rank",    1,      dimid_rank  ) )
      
      dimid_2d = (/ dimid_corner, dimid_column /) 

      ! Define the output variables
      call check( nf90_def_var(ncid, "grid_area",       NF90_DOUBLE, dimid_column, varid(1) ) )
      call check( nf90_def_var(ncid, "grid_center_lat", NF90_DOUBLE, dimid_column, varid(2) ) )
      call check( nf90_def_var(ncid, "grid_center_lon", NF90_DOUBLE, dimid_column, varid(3) ) )
      call check( nf90_def_var(ncid, "grid_corner_lat", NF90_DOUBLE, dimid_2d,     varid(4) ) )
      call check( nf90_def_var(ncid, "grid_corner_lon", NF90_DOUBLE, dimid_2d,     varid(5) ) )
      call check( nf90_def_var(ncid, "grid_imask",      NF90_INT,    dimid_column, varid(6) ) )
      call check( nf90_def_var(ncid, "grid_dims",       NF90_INT,    dimid_rank,   varid(7) ) )

      ! Add units attribute
      call check( nf90_put_att(ncid, varid(1), "units", "radians^2") )
      call check( nf90_put_att(ncid, varid(2), "units", "degrees") )
      call check( nf90_put_att(ncid, varid(3), "units", "degrees") )
      call check( nf90_put_att(ncid, varid(4), "units", "degrees") )
      call check( nf90_put_att(ncid, varid(5), "units", "degrees") )

      ! End definition mode
      call check( nf90_enddef(ncid) )

      ! Write the actual data to the file   
      call check( nf90_put_var(ncid, varid(1), area(:) ) )
      call check( nf90_put_var(ncid, varid(2), center_lat_deg(:)   ) )
      call check( nf90_put_var(ncid, varid(3), center_lon_deg(:)   ) )
      call check( nf90_put_var(ncid, varid(4), corner_lat_deg_alt(:,:) ) )
      call check( nf90_put_var(ncid, varid(5), corner_lon_deg_alt(:,:) ) )
      call check( nf90_put_var(ncid, varid(6), mask(:) ) )
      call check( nf90_put_var(ncid, varid(7), 1 ) )

      ! Close the file
      call check( nf90_close(ncid) )

   end if ! masterproc
   !----------------------------------------------------------------------------

   ! call endrun("write_grid_scrip: Exiting after writing scrip file: "//trim(filename) )

end subroutine write_grid_scrip 
!===================================================================================================
!===================================================================================================
end module write_grid_scrip_mod