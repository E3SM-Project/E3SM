MODULE mapread_mod

   use maptype_mod
   use shr_sys_mod

   implicit none

#include <netcdf.inc>

!===============================================================================
CONTAINS
!===============================================================================

SUBROUTINE mapread_dest_grid(map, filename)

   !--- modules ---

   implicit none

   !--- includes ---

   !--- arguments ---
   type(sMatrix), intent(inout) :: map       ! sMatrix info to be read in
   character(*) , intent(in)    :: filename  ! name of data file

   !--- local ---
   character(strLen)     :: units     ! netCDF attribute name string
   integer               :: rcode   ! netCDF routine return code
   integer               :: fid     ! netCDF file      ID
   integer               :: vid     ! netCDF variable  ID
   integer               :: did     ! netCDF dimension ID
   integer               :: grid_rank
   integer, allocatable, dimension(:) :: grid_dims

   character(*),parameter :: subName = "(mapread_dst_grid) "
   !--- formats ---
   character(len=*),parameter :: F00 = "('(mapread_dst_grid) ',3a)"
   character(len=*),parameter :: F02 = "('(mapread_dst_grid) ',a11,a3,60(a1))"

   write(6,F00) 'ocn data file',' = ',trim(filename)
   rcode = nf_open(filename,NF_NOWRITE,fid)

   rcode = nf_inq_dimid (fid, 'grid_size' , did)
   rcode = nf_inq_dimlen(fid, did   , map%n_b  )
   rcode = nf_inq_dimid (fid, 'grid_corners' , did)
   rcode = nf_inq_dimlen(fid, did   , map%nv_b  )
   rcode = nf_inq_dimid (fid, 'grid_rank', did)
   rcode = nf_inq_dimlen(fid, did   , grid_rank)
   allocate(grid_dims(grid_rank))
   rcode = nf_inq_varid  (fid, 'grid_dims', vid)
   rcode = nf_get_var_int(fid, vid   , grid_dims)
   if (grid_rank.eq.1) then
     map%ni_b = grid_dims(1)
     map%nj_b = 1
   elseif (grid_rank.eq.2) then
     map%ni_b = grid_dims(1)
     map%nj_b = grid_dims(2)
   else
      deallocate(grid_dims)
      write(6,*) 'ERROR: grid_rank is ',grid_rank,' in ',trim(filename)
      call shr_sys_abort(subName//"ERROR: filename grid_rank")
   endif
   deallocate(grid_dims)
   map%dims_b(1) = map%ni_b
   map%dims_b(2) = map%nj_b

   allocate(map%  xc_b(         map%n_b)) ! x-coordinates of center
   allocate(map%  yc_b(         map%n_b)) ! y-coordinates of center
   allocate(map%  xv_b(map%nv_b,map%n_b)) ! x-coordinates of verticies
   allocate(map%  yv_b(map%nv_b,map%n_b)) ! y-coordinates of verticies
   allocate(map%mask_b(         map%n_b)) ! domain mask
   allocate(map%area_b(         map%n_b)) ! grid cell area
   allocate(map%frac_b(         map%n_b)) ! grid cell area
   allocate(map%sn1            (map%n_b))
   allocate(map%sn2            (map%n_b))

   rcode = nf_inq_varid     (fid,'grid_center_lon'  ,vid)
   rcode = nf_get_var_double(fid,vid     ,map%xc_b )
   units = "" ! units needs to be emptied before reading from netCDF file
   rcode = nf_get_att_text(fid, vid, "units", units)
   if (trim(units).eq."radians") then
      map%xc_b = map%xc_b * RADtoDEG
   end if

   rcode = nf_inq_varid     (fid,'grid_center_lat'  ,vid)
   rcode = nf_get_var_double(fid,vid     ,map%yc_b )
   units = "" ! units needs to be emptied before reading from netCDF file
   rcode = nf_get_att_text(fid, vid, "units", units)
   if (trim(units).eq."radians") then
      map%yc_b = map%yc_b * RADtoDEG
   end if

   rcode = nf_inq_varid     (fid,'grid_corner_lon'  ,vid)
   rcode = nf_get_var_double(fid,vid     ,map%xv_b )
   units = "" ! units needs to be emptied before reading from netCDF file
   rcode = nf_get_att_text(fid, vid, "units", units)
   if (trim(units).eq."radians") then
      map%xv_b = map%xv_b * RADtoDEG
   end if

   rcode = nf_inq_varid     (fid,'grid_corner_lat'  ,vid)
   rcode = nf_get_var_double(fid,vid     ,map%yv_b )
   units = "" ! units needs to be emptied before reading from netCDF file
   rcode = nf_get_att_text(fid, vid, "units", units)
   if (trim(units).eq."radians") then
      map%yv_b = map%yv_b * RADtoDEG
   end if

   rcode = nf_inq_varid     (fid,'grid_imask',vid )
   rcode = nf_get_var_int   (fid,vid     ,map%mask_b)
   rcode = nf_inq_varid     (fid,'grid_area',vid )
   if (rcode.eq.0) then
      rcode = nf_get_var_double(fid,vid     ,map%area_b)
   else
      write(6,*) "ERROR: could not find variable grid_area in destination grid input file!"
      stop
   end if

   map%frac_b = map%mask_b * 1.0_r8
   map%domain_b = trim(filename)
   rcode = nf_close(fid)

END SUBROUTINE mapread_dest_grid

!===============================================================================
!===============================================================================
END MODULE mapread_mod
!===============================================================================
