!===============================================================================
! SVN $Id: map_mod.F90 46983 2013-05-09 22:08:12Z tcraig $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/tools/mapping/trunk_tags/mapping_130509/gen_mapping_files/runoff_to_ocn/src/map_mod.F90 $
!===============================================================================

MODULE map_mod

   use shr_sys_mod
   use shr_timer_mod
   use kind_mod

   implicit none

#include <netcdf.inc>

   integer,parameter :: strLen = 240
   integer,parameter :: nv     = 4

   real(r8)   ,parameter :: pi =  3.14159265358979323846_r8
   real(r8)   ,parameter :: rEarth     =  6.37122e+6         ! radius of earth (m)
   real(r8)   ,parameter :: DEGtoRAD = pi/180.0_r8  ! degrees to radians
   real(r8)   ,parameter :: RADtoDEG = 180.0_r8/pi  ! radians to degrees

   integer,parameter :: nibx = 360          ! size of sort bin boxes
   integer,parameter :: njbx = nibx/2
   integer,parameter :: vcells_req = 20     ! number of valid cells required per bin

   !----------------------------------------------------------------------------
   ! sparse matrix data type
   !----------------------------------------------------------------------------
   TYPE sMatrix

     !--- global text attributes ---
     character(strLen) :: title
     character(strLen) :: normal
     character(strLen) :: method
     character(strLen) :: history
     character(strLen) :: convention
     character(strLen) :: domain_a
     character(strLen) :: domain_b

     !--- domain a ---
     integer         ::    n_a      ! number of non-zero matrix elements
     integer         ::   ni_a      ! number of 2d array i indicies
     integer         ::   nj_a      ! number of 2d array j indicies
     real(r8)   ,pointer ::   xc_a(:)   ! x-coords of centers   ~ deg east
     real(r8)   ,pointer ::   yc_a(:)   ! y-coords of centers   ~ deg north
     real(r8)   ,pointer ::   xv_a(:,:) ! x-coords of verticies ~ deg east, (nv,n)
     real(r8)   ,pointer ::   yv_a(:,:) ! y-coords of verticies ~ deg north (nv,n)
     integer,pointer :: mask_a(:)   ! mask: 0 <=> out-of-domain (invalid data)
     real(r8)   ,pointer :: area_a(:)   ! area of grid cell ~ radians squared
     integer         :: dims_a(2)       ! hardwire to 2 for now

     !--- domain b ---
     integer         ::    n_b      ! number of non-zero matrix elements
     integer         ::   ni_b      ! number of 2d array i indicies
     integer         ::   nj_b      ! number of 2d array j indicies
     real(r8)   ,pointer ::   xc_b(:)   ! x-coords of centers   ~ deg east
     real(r8)   ,pointer ::   yc_b(:)   ! y-coords of centers   ~ deg north
     real(r8)   ,pointer ::   xv_b(:,:) ! x-coords of verticies ~ deg east, (nv,n)
     real(r8)   ,pointer ::   yv_b(:,:) ! y-coords of verticies ~ deg north (nv,n)
     integer,pointer :: mask_b(:)   ! mask: 0 <=> out-of-domain (invalid data)
     real(r8)   ,pointer :: area_b(:)   ! area of grid cell ~ radians squared
     integer         :: dims_b(2)       ! hardwire to 2 for now

     !--- fraction of cell mapped to domain b or from domain a ---
     real(r8)   ,pointer :: frac_a(:)   ! area of grid cell ~ radians squared
     real(r8)   ,pointer :: frac_b(:)   ! area of grid cell ~ radians squared

     !--- map: a->b ---
     integer         :: n_s         ! number of non-zero matrix elements
     real(kind=r8)   ,pointer :: s  (:)      ! the non-zero matrix elements
     integer,pointer :: row(:)      ! matrix row corresponding to each element
     integer,pointer :: col(:)      ! matrix col corresponding to each element

     !--- used for OMP/threading in mat-mult ---
     integer,pointer :: sn1(:)      ! # links in a given row
     integer,pointer :: sn2(:)      ! # links previous to a given row

     !--- required for computing NN maps ---
     integer :: imin_b(nibx,njbx)   ! xc_b least index for lat 0:360
     integer :: imax_b(nibx,njbx)   ! xc_b max index for lat 0:360
     integer :: jmin_b(nibx,njbx)   ! yc_b least index for lat 0:90
     integer :: jmax_b(nibx,njbx)   ! yc_b max index for lat 0:90

   END TYPE sMatrix

   SAVE

!===============================================================================
CONTAINS 
!===============================================================================

SUBROUTINE map_read(map, filename)

   !--- modules ---

   implicit none

   !--- includes ---

   !--- arguments ---
   type(sMatrix), intent(inout) :: map       ! sMatrix info to be read in
   character(*) , intent(in)    :: filename  ! name of data file

   !--- local ---
   integer         :: i,n,m         ! generic indicies

   character(strLen)     :: str     ! variable length char string
   character(strLen)     :: attstr  ! netCDF attribute name string
   integer               :: rcode   ! netCDF routine return code
   integer               :: fid     ! netCDF file      ID
   integer               :: vid     ! netCDF variable  ID
   integer               :: did     ! netCDF dimension ID

   !--- formats ---
   character(len=*),parameter :: F00 = "('(map_read) ',3a)"
   character(len=*),parameter :: F02 = "('(map_read) ',a11,a3,60(a1))"

!-------------------------------------------------------------------------------
! PURPOSE:
! o reads map matrix information from netCDF data file
! 
! NOTE:
!-------------------------------------------------------------------------------

   !-----------------------------------------------
   ! get global attributes
   !-----------------------------------------------
   write(6,F00) 'data file',' = ',trim(filename)
   rcode = nf_open(filename,NF_NOWRITE,fid)
   if (rcode.ne.NF_NOERR) then
      write(6,F00) ' ERROR in filename ',trim(filename)
      call shr_sys_abort()
   endif
   do i=1,7
     if (i == 1) attstr = 'title'      
     if (i == 2) attstr = 'normalization'
     if (i == 3) attstr = 'map_method'
     if (i == 4) attstr = 'conventions'
     if (i == 5) attstr = 'history'
     if (i == 6) attstr = 'domain_a'
     if (i == 7) attstr = 'domain_b'

     rcode = nf_inq_attlen  (fid, NF_GLOBAL, trim(attstr), n)
     rcode = nf_get_att_text(fid, NF_GLOBAL, trim(attstr), str)
     write(6,F00) trim(attstr),' = ',str(1:n)

     if (i == 1) map%title      = str(1:n)
     if (i == 2) map%normal     = str(1:n)
     if (i == 3) map%method     = str(1:n)
     if (i == 4) map%convention = str(1:n)
     if (i == 5) map%history    = str(1:n)
     if (i == 6) map%domain_a   = str(1:n)
     if (i == 7) map%domain_b   = str(1:n)

   end do

   !-----------------------------------------------
   ! get "a" domain info 
   !-----------------------------------------------
   rcode = nf_inq_dimid (fid, 'n_a' , did)
   rcode = nf_inq_dimlen(fid, did   , map%n_a  )
   map%dims_a(1) = map%n_a
   map%dims_a(2) = 1
   rcode = nf_inq_dimid (fid, 'ni_a', did)
   if (rcode.eq.NF_NOERR) then
      rcode = nf_inq_dimlen(fid, did   , map%ni_a )
      map%dims_a(1) = map%ni_a
   endif
   rcode = nf_inq_dimid (fid, 'nj_a', did)
   if (rcode.eq.NF_NOERR) then
      rcode = nf_inq_dimlen(fid, did   , map%nj_a )
      map%dims_a(2) = map%nj_a
   endif

   allocate(map%  xc_a(   map%n_a)) ! x-coordinates of center
   allocate(map%  yc_a(   map%n_a)) ! y-coordinates of center
   allocate(map%  xv_a(nv,map%n_a)) ! x-coordinates of verticies
   allocate(map%  yv_a(nv,map%n_a)) ! y-coordinates of verticies
   allocate(map%mask_a(   map%n_a)) ! domain mask
   allocate(map%area_a(   map%n_a)) ! grid cell area
   allocate(map%frac_a(   map%n_a)) ! grid cell area

   rcode = nf_inq_varid     (fid,'src_grid_dims',vid )
   if (rcode.eq.NF_NOERR) then
      rcode = nf_get_var_int   (fid,vid     ,map%dims_a)
   endif
   rcode = nf_inq_varid     (fid,'xc_a'  ,vid)
   rcode = nf_get_var_double(fid,vid     ,map%xc_a )
   rcode = nf_inq_varid     (fid,'yc_a'  ,vid)
   rcode = nf_get_var_double(fid,vid     ,map%yc_a )
   rcode = nf_inq_varid     (fid,'xv_a'  ,vid)
   rcode = nf_get_var_double(fid,vid     ,map%xv_a )
   rcode = nf_inq_varid     (fid,'yv_a'  ,vid)
   rcode = nf_get_var_double(fid,vid     ,map%yv_a )
   rcode = nf_inq_varid     (fid,'mask_a',vid )
   rcode = nf_get_var_int   (fid,vid     ,map%mask_a)
   rcode = nf_inq_varid     (fid,'area_a',vid )
   rcode = nf_get_var_double(fid,vid     ,map%area_a)
   rcode = nf_inq_varid     (fid,'frac_a',vid )
   rcode = nf_get_var_double(fid,vid     ,map%frac_a)

   !-----------------------------------------------
   ! get "b" domain info 
   !-----------------------------------------------
   rcode = nf_inq_dimid (fid, 'n_b' , did)
   rcode = nf_inq_dimlen(fid, did   , map%n_b  )
   map%dims_b(1) = map%n_b
   map%dims_b(2) = 1
   rcode = nf_inq_dimid (fid, 'ni_b', did)
   if (rcode.eq.NF_NOERR) then
      rcode = nf_inq_dimlen(fid, did   , map%ni_b )
      map%dims_b(1) = map%ni_b
   endif
   rcode = nf_inq_dimid (fid, 'nj_b', did)
   if (rcode.eq.NF_NOERR) then
      rcode = nf_inq_dimlen(fid, did   , map%nj_b )
      map%dims_b(2) = map%nj_b
   endif

   allocate(map%  xc_b(   map%n_b)) ! x-coordinates of center
   allocate(map%  yc_b(   map%n_b)) ! y-coordinates of center
   allocate(map%  xv_b(nv,map%n_b)) ! x-coordinates of verticies
   allocate(map%  yv_b(nv,map%n_b)) ! y-coordinates of verticies
   allocate(map%mask_b(   map%n_b)) ! domain mask
   allocate(map%area_b(   map%n_b)) ! grid cell area
   allocate(map%frac_b(   map%n_b)) ! grid cell area

   rcode = nf_inq_varid     (fid,'dst_grid_dims',vid )
   if (rcode.eq.NF_NOERR) then
      rcode = nf_get_var_int   (fid,vid     ,map%dims_b)
   endif
   rcode = nf_inq_varid     (fid,'xc_b'  ,vid)
   rcode = nf_get_var_double(fid,vid     ,map%xc_b )
   rcode = nf_inq_varid     (fid,'yc_b'  ,vid)
   rcode = nf_get_var_double(fid,vid     ,map%yc_b )
   rcode = nf_inq_varid     (fid,'xv_b'  ,vid)
   rcode = nf_get_var_double(fid,vid     ,map%xv_b )
   rcode = nf_inq_varid     (fid,'yv_b'  ,vid)
   rcode = nf_get_var_double(fid,vid     ,map%yv_b )
   rcode = nf_inq_varid     (fid,'mask_b',vid )
   rcode = nf_get_var_int   (fid,vid     ,map%mask_b)
   rcode = nf_inq_varid     (fid,'area_b',vid )
   rcode = nf_get_var_double(fid,vid     ,map%area_b)
   rcode = nf_inq_varid     (fid,'frac_b',vid )
   rcode = nf_get_var_double(fid,vid     ,map%frac_b)

   !-----------------------------------------------
   ! get matrix info 
   !-----------------------------------------------
   rcode = nf_inq_dimid (fid, 'n_s', did)  ! size of sparse matrix
   rcode = nf_inq_dimlen(fid, did  , map%n_s)
   allocate(map%s  (map%n_s))
   allocate(map%row(map%n_s))
   allocate(map%col(map%n_s))
   allocate(map%sn1(map%n_b))
   allocate(map%sn2(map%n_b))
   rcode = nf_inq_varid     (fid,'S'  ,vid)
   rcode = nf_get_var_double(fid,vid  ,map%s  )
   rcode = nf_inq_varid     (fid,'row',vid)
   rcode = nf_get_var_int   (fid,vid  ,map%row)
   rcode = nf_inq_varid     (fid,'col',vid)
   rcode = nf_get_var_int   (fid,vid  ,map%col)

   rcode = nf_close(fid)

END SUBROUTINE map_read

!===============================================================================
!===============================================================================

SUBROUTINE map_gridRead(map, rfilename, ofilename, gridtype)

   !--- modules ---

   implicit none

   !--- includes ---

   !--- arguments ---
   type(sMatrix), intent(inout) :: map       ! sMatrix info to be read in
   character(*) , intent(in)    :: rfilename ! name of rtm rdirc file
   character(*) , intent(in)    :: ofilename ! name of ocn scrip grid file
   character(*) , intent(in)    :: gridtype  ! type of roff grid data 

   !--- local ---
   integer         :: i,j,n,m         ! generic indicies

   character(strLen)     :: str     ! variable length char string
   character(strLen)     :: attstr  ! netCDF attribute name string
   integer               :: rcode   ! netCDF routine return code
   integer               :: fid     ! netCDF file      ID
   integer               :: vid     ! netCDF variable  ID
   integer               :: did     ! netCDF dimension ID
   integer               :: grid_rank
   integer               :: grid_dims(2)
   real(r8)              :: lat,lon,rdirc
   integer               :: ilat,ilon
   integer               :: nm,np,cnt,ibx,jbx,ibxj,jbxj,ni,nj,nij
   integer               :: im,ip,jp,jmax,jmin,imin,imax,jm,ibxi,jbxi,ibxij,jbxij
   integer               :: ispval,ihold,i1,j1,i2,j2,cnt2,i3,j3
   integer               :: vcells
   integer,allocatable   :: timin(:,:),timax(:,:),tjmin(:,:),tjmax(:,:)
   logical               :: lxv_a, lyv_a, lmask_a, larea_a

   !--- formats ---
   character(*),parameter :: subName = "(map_gridRead) "
   character(*),parameter :: F00 = "('(map_gridRead) ',3a)"
   character(*),parameter :: F01 = "('(map_gridRead) ',a,2i8)"
   character(*),parameter :: F02 = "('(map_gridRead) ',a11,a3,60(a1))"

!-------------------------------------------------------------------------------
! PURPOSE:
! o reads map matrix information from netCDF data file
! 
! NOTE:
!-------------------------------------------------------------------------------

   !-----------------------------------------------
   ! set global attributesa
   !-----------------------------------------------

   map%title      = 'gen_runoff mapping generation'
   map%normal     = 'conservative'
   map%method     = 'nearest neighbor smoothing'
   map%convention = 'NCAR-CCSM'
   map%history    = 'history'
   map%domain_a   = trim(rfilename)
   map%domain_b   = trim(ofilename)

   if (trim(gridtype) == "rtm") then
      !-------------------------------------------------------------------------
      write(*,F00) "read source domain info -- hardwired for 1/2 degree clm/rtm rdirc ascii file"
      !-------------------------------------------------------------------------

      map%ni_a = 720
      map%nj_a = 360
      map%n_a  = map%ni_a * map%nj_a
      map%dims_a(1) = map%ni_a
      map%dims_a(2) = map%nj_a

      allocate(map%  xc_a(   map%n_a)) ! x-coordinates of center
      allocate(map%  yc_a(   map%n_a)) ! y-coordinates of center
      allocate(map%  xv_a(nv,map%n_a)) ! x-coordinates of verticies
      allocate(map%  yv_a(nv,map%n_a)) ! y-coordinates of verticies
      allocate(map%mask_a(   map%n_a)) ! domain mask
      allocate(map%area_a(   map%n_a)) ! grid cell area
      allocate(map%frac_a(   map%n_a)) ! grid cell area

      fid = 22
      open(fid,file=trim(rfilename))

      do n = 1,map%n_a
         read(fid,*) lat,lon,rdirc
         map%xc_a(  n) = lon
         map%yc_a(  n) = lat

         map%xv_a(1,n) = lon - 0.25
         map%xv_a(2,n) = lon + 0.25
         map%xv_a(3,n) = lon + 0.25
         map%xv_a(4,n) = lon - 0.25

         map%yv_a(1,n) = lat - 0.25
         map%yv_a(2,n) = lat - 0.25
         map%yv_a(3,n) = lat + 0.25
         map%yv_a(4,n) = lat + 0.25

         map%area_a(n) = 0.5 * 0.5 * cos(lat*DEGtoRAD) * DEGtoRAD * DEGtoRAD

         if (abs(rdirc) < 0.5) then
            map%mask_a(n) = 1
            map%frac_a(n) = 1.0_r8
         else
            map%mask_a(n) = 0
            map%frac_a(n) = 0.0_r8
         endif
      enddo

      close(fid)
   else if (trim(gridtype) == "obs") then
      !-------------------------------------------------------------------------
      write(*,F00) "read source domain info -- Dai/Trenberth observed runoff"
      !-------------------------------------------------------------------------
      write(6,F00) 'runoff data file',' = ',trim(rfilename)
      rcode = nf_open(rfilename,NF_NOWRITE,fid)

      rcode = nf_inq_dimid (fid, 'ni' , did)
      if (rcode.ne.NF_NOERR) then
        write(6,F00) ' ERROR in filename ',trim(rfilename),': can not find dim ni'
        call shr_sys_abort()
      endif
      rcode = nf_inq_dimlen(fid, did  , ni )

      rcode = nf_inq_dimid (fid, 'nj' , did)
      if (rcode.ne.NF_NOERR) then
        write(6,F00) ' ERROR in filename ',trim(rfilename),': can not find dim nj'
        call shr_sys_abort()
      endif
      rcode = nf_inq_dimlen(fid, did  , nj )

      write(6,F01) 'dimensions ni,nj = ',ni,nj

      map%ni_a = ni
      map%nj_a = nj
      map%n_a  = ni*nj
      map%dims_a(1) = map%ni_a
      map%dims_a(2) = map%nj_a

      allocate(map%  xc_a(   map%n_a)) ! x-coordinates of center
      allocate(map%  yc_a(   map%n_a)) ! y-coordinates of center
      allocate(map%  xv_a(nv,map%n_a)) ! x-coordinates of verticies
      allocate(map%  yv_a(nv,map%n_a)) ! y-coordinates of verticies
      allocate(map%mask_a(   map%n_a)) ! domain mask
      allocate(map%area_a(   map%n_a)) ! grid cell area
      allocate(map%frac_a(   map%n_a)) ! frac of cell mapped to dest

      rcode = nf_inq_varid     (fid,'xc'  ,vid)
      if (rcode.ne.NF_NOERR) then
        write(6,*) "ERROR: could not find variable xc in input file!"
        call shr_sys_abort()
      endif
      rcode = nf_get_var_double(fid,vid   ,map%xc_a )

      rcode = nf_inq_varid     (fid,'yc'  ,vid)
      if (rcode.ne.NF_NOERR) then
        write(6,*) "ERROR: could not find variable yc in input file!"
        call shr_sys_abort()
      endif
      rcode = nf_get_var_double(fid,vid   ,map%yc_a )

      rcode = nf_inq_varid     (fid,'xv'  ,vid)
      if (rcode.eq.0) then
        lxv_a = .true.
        rcode = nf_get_var_double(fid,vid   ,map%xv_a )
      else
        lxv_a = .false.
        write(6,*) "WARNING: could not find variable xv in input file!"
      endif

      rcode = nf_inq_varid     (fid,'yv'  ,vid)
      if (rcode.eq.0) then
        lyv_a = .true.
        rcode = nf_get_var_double(fid,vid   ,map%yv_a )
      else
        lyv_a = .false.
        write(6,*) "WARNING: could not find variable yv in input file!"
      endif

      rcode = nf_inq_varid     (fid,'mask',vid )
      if (rcode.eq.0) then
        lmask_a = .true.
        rcode = nf_get_var_int   (fid,vid   ,map%mask_a)
      else
        lmask_a = .false.
        write(6,*) "WARNING: could not find variable mask in input file!"
      endif

      rcode = nf_inq_varid     (fid,'area',vid )
      if (rcode.eq.0) then
        larea_a = .true.
        rcode = nf_get_var_double(fid,vid   ,map%area_a)
      else
        larea_a = .false.
        write(6,*) "WARNING: could not find variable area in input file!"
      endif


      rcode = nf_close(fid)

      do n = 1,map%n_a
         lon = map%xc_a(n)
         lat = map%yc_a(n)

         if (.not.lxv_a) then
           ! Hardcoded for r0.5
           map%xv_a(1,n) = lon - 0.25_r8
           map%xv_a(2,n) = lon + 0.25_r8
           map%xv_a(3,n) = lon + 0.25_r8
           map%xv_a(4,n) = lon - 0.25_r8
         end if

         if (.not.lyv_a) then
           ! Hardcoded for r0.5
           map%yv_a(1,n) = lat - 0.25_r8
           map%yv_a(2,n) = lat - 0.25_r8
           map%yv_a(3,n) = lat + 0.25_r8
           map%yv_a(4,n) = lat + 0.25_r8
         end if

         if (.not.larea_a) then
           ! Hardcoded for r0.5
           map%area_a(n) = 0.5*0.5*cos(lat*DEGtoRAD)*DEGtoRAD*DEGtoRAD
         end if

         if (.not.lmask_a) then
           map%mask_a(n) = 1
         end if

         if (map%mask_a(n) .ne. 0) then
            map%mask_a(n) = 1
            map%frac_a(n) = 1.0_r8
         else
            map%mask_a(n) = 0
            map%frac_a(n) = 0.0_r8
         endif
      enddo

!     !--- safe to assume units are already degrees? ---
!     if (units == radians)
!     map%xc_b = map%xc_b * RADtoDEG 
!     map%yc_b = map%yc_b * RADtoDEG
!     map%xv_b = map%xv_b * RADtoDEG
!     map%yv_b = map%yv_b * RADtoDEG

   else if (trim(gridtype) == "scrip") then

      !----------------------------------------------------------------------------
      write(*,F00) "read source domain info -- scrip grid"
      !----------------------------------------------------------------------------
      write(6,F00) 'rof data file',' = ',trim(rfilename)
      rcode = nf_open(rfilename,NF_NOWRITE,fid)

      rcode = nf_inq_dimid (fid, 'grid_size' , did)
      rcode = nf_inq_dimlen(fid, did   , map%n_a  )
      rcode = nf_inq_dimid (fid, 'grid_rank', did)
      rcode = nf_inq_dimlen(fid, did   , grid_rank)
      if (grid_rank /= 2) then
         write(6,*) 'ERROR: grid_rank is ',grid_rank,' in ',trim(rfilename)
         call shr_sys_abort(subName//"ERROR: rfilename grid_rank")
      endif
      rcode = nf_inq_varid  (fid, 'grid_dims', vid)
      rcode = nf_get_var_int(fid, vid   , grid_dims)
      map%ni_a = grid_dims(1)
      map%nj_a = grid_dims(2)
      map%dims_a(1) = map%ni_a
      map%dims_a(2) = map%nj_a

      allocate(map%  xc_a(   map%n_a)) ! x-coordinates of center
      allocate(map%  yc_a(   map%n_a)) ! y-coordinates of center
      allocate(map%  xv_a(nv,map%n_a)) ! x-coordinates of verticies
      allocate(map%  yv_a(nv,map%n_a)) ! y-coordinates of verticies
      allocate(map%mask_a(   map%n_a)) ! domain mask
      allocate(map%area_a(   map%n_a)) ! grid cell area
      allocate(map%frac_a(   map%n_a)) ! grid cell area

      rcode = nf_inq_varid     (fid,'grid_center_lon'  ,vid)
      rcode = nf_get_var_double(fid,vid     ,map%xc_a )
      rcode = nf_inq_varid     (fid,'grid_center_lat'  ,vid)
      rcode = nf_get_var_double(fid,vid     ,map%yc_a )
      rcode = nf_inq_varid     (fid,'grid_corner_lon'  ,vid)
      rcode = nf_get_var_double(fid,vid     ,map%xv_a )
      rcode = nf_inq_varid     (fid,'grid_corner_lat'  ,vid)
      rcode = nf_get_var_double(fid,vid     ,map%yv_a )
      rcode = nf_inq_varid     (fid,'grid_imask',vid )
      rcode = nf_get_var_int   (fid,vid     ,map%mask_a)
      rcode = nf_inq_varid     (fid,'grid_area',vid )
      rcode = nf_get_var_double(fid,vid     ,map%area_a)

      map%xc_a = map%xc_a * RADtoDEG
      map%yc_a = map%yc_a * RADtoDEG
      map%xv_a = map%xv_a * RADtoDEG
      map%yv_a = map%yv_a * RADtoDEG
      map%frac_a = map%mask_a * 1.0_r8

      rcode = nf_close(fid)

   else
      !-------------------------------------------------------------------------
      ! no other valid choices
      !-------------------------------------------------------------------------
      write(6,F00) 'ERROR: no implementation for gridtype = ',trim(gridtype)
      stop
   end if
 
   !----------------------------------------------------------------------------
   write(*,F00) "read destination domain info -- pop grid"
   !----------------------------------------------------------------------------
   write(6,F00) 'ocn data file',' = ',trim(ofilename)
   rcode = nf_open(ofilename,NF_NOWRITE,fid)

   rcode = nf_inq_dimid (fid, 'grid_size' , did)
   rcode = nf_inq_dimlen(fid, did   , map%n_b  )
   rcode = nf_inq_dimid (fid, 'grid_rank', did)
   rcode = nf_inq_dimlen(fid, did   , grid_rank)
   if (grid_rank /= 2) then
      write(6,*) 'ERROR: grid_rank is ',grid_rank,' in ',trim(ofilename)
      call shr_sys_abort(subName//"ERROR: ofilename grid_rank")
   endif
   rcode = nf_inq_varid  (fid, 'grid_dims', vid)
   rcode = nf_get_var_int(fid, vid   , grid_dims)
   map%ni_b = grid_dims(1)
   map%nj_b = grid_dims(2)
   map%dims_b(1) = map%ni_b
   map%dims_b(2) = map%nj_b

   allocate(map%  xc_b(   map%n_b)) ! x-coordinates of center
   allocate(map%  yc_b(   map%n_b)) ! y-coordinates of center
   allocate(map%  xv_b(nv,map%n_b)) ! x-coordinates of verticies
   allocate(map%  yv_b(nv,map%n_b)) ! y-coordinates of verticies
   allocate(map%mask_b(   map%n_b)) ! domain mask
   allocate(map%area_b(   map%n_b)) ! grid cell area
   allocate(map%frac_b(   map%n_b)) ! grid cell area
   allocate(map%sn1      (map%n_b))
   allocate(map%sn2      (map%n_b))

   rcode = nf_inq_varid     (fid,'grid_center_lon'  ,vid)
   rcode = nf_get_var_double(fid,vid     ,map%xc_b )
   rcode = nf_inq_varid     (fid,'grid_center_lat'  ,vid)
   rcode = nf_get_var_double(fid,vid     ,map%yc_b )
   rcode = nf_inq_varid     (fid,'grid_corner_lon'  ,vid)
   rcode = nf_get_var_double(fid,vid     ,map%xv_b )
   rcode = nf_inq_varid     (fid,'grid_corner_lat'  ,vid)
   rcode = nf_get_var_double(fid,vid     ,map%yv_b )
   rcode = nf_inq_varid     (fid,'grid_imask',vid )
   rcode = nf_get_var_int   (fid,vid     ,map%mask_b)
   rcode = nf_inq_varid     (fid,'grid_area',vid )
   rcode = nf_get_var_double(fid,vid     ,map%area_b)

   map%xc_b = map%xc_b * RADtoDEG
   map%yc_b = map%yc_b * RADtoDEG
   map%xv_b = map%xv_b * RADtoDEG
   map%yv_b = map%yv_b * RADtoDEG
   map%frac_b = map%mask_b * 1.0_r8

   rcode = nf_close(fid)

   !----------------------------------------------------------------------------
   write(*,F00) "derive info required to compute NN map"
   !----------------------------------------------------------------------------

   ispval = -999
   map%imin_b = 5*map%ni_b
   map%imax_b = ispval
   map%jmin_b = 5*map%nj_b
   map%jmax_b = ispval

   do j = 1,map%nj_b-1
   do i = 1,map%ni_b
      jm = j
      jp = j+1
      im = i
      ip = mod(i,map%ni_b) + 1

      n   = (jm-1)*map%ni_b + im
      nj  = (jp-1)*map%ni_b + im
      ni  = (jm-1)*map%ni_b + ip
      nij = (jp-1)*map%ni_b + ip

      call map_bxindex(map%xc_b(n  ),map%yc_b(n  ),ibx  ,jbx  )
      call map_bxindex(map%xc_b(nj ),map%yc_b(nj ),ibxj ,jbxj )
      call map_bxindex(map%xc_b(ni ),map%yc_b(ni ),ibxi ,jbxi )
      call map_bxindex(map%xc_b(nij),map%yc_b(nij),ibxij,jbxij)

      imin = min(ibx,ibxi,ibxj,ibxij)
      imax = max(ibx,ibxi,ibxj,ibxij)
      jmin = min(jbx,jbxi,jbxj,jbxij)
      jmax = max(jbx,jbxi,jbxj,jbxij)

      ! --- if wraparound in bx, swap index order
      if (imax - imin > nibx/2) then
         ihold = imin
         imin  = imax
         imax  = ihold
      endif

      ip = mod(ip-1,map%ni_b) + 1

      if (imax < imin) then
         do j1 = jmin,jmax
         do i1 = 1,imax
            map%jmin_b(i1,j1) =  min(map%jmin_b(i1,j1),jm)
            map%jmax_b(i1,j1) =  max(map%jmax_b(i1,j1),jp)
            if (map%imax_b(i1,j1) == ispval) then
               map%imin_b(i1,j1) = im
               map%imax_b(i1,j1) = ip
            else
               if (abs(map%imin_b(i1,j1) - im) > nibx/2) then
                  map%imin_b(i1,j1) = max(map%imin_b(i1,j1),im)
               else
                  map%imin_b(i1,j1) = min(map%imin_b(i1,j1),im)
               endif
               if (abs(map%imax_b(i1,j1) - ip) > nibx/2) then
                  map%imax_b(i1,j1) = min(map%imax_b(i1,j1),ip)
               else
                  map%imax_b(i1,j1) = max(map%imax_b(i1,j1),ip)
               endif
            endif
         enddo
         do i1 = imin,nibx
            map%jmin_b(i1,j1) =  min(map%jmin_b(i1,j1),jm)
            map%jmax_b(i1,j1) =  max(map%jmax_b(i1,j1),jp)
            if (map%imax_b(i1,j1) == ispval) then
               map%imin_b(i1,j1) = im
               map%imax_b(i1,j1) = ip
            else
               if (abs(map%imin_b(i1,j1) - im) > nibx/2) then
                  map%imin_b(i1,j1) = max(map%imin_b(i1,j1),im)
               else
                  map%imin_b(i1,j1) = min(map%imin_b(i1,j1),im)
               endif
               if (abs(map%imax_b(i1,j1) - ip) > nibx/2) then
                  map%imax_b(i1,j1) = min(map%imax_b(i1,j1),ip)
               else
                  map%imax_b(i1,j1) = max(map%imax_b(i1,j1),ip)
               endif
            endif
         enddo
         enddo
      else
         do j1 = jmin,jmax
         do i1 = imin,imax
            map%jmin_b(i1,j1) =  min(map%jmin_b(i1,j1),jm)
            map%jmax_b(i1,j1) =  max(map%jmax_b(i1,j1),jp)
            if (map%imax_b(i1,j1) == ispval) then
               map%imin_b(i1,j1) = im
               map%imax_b(i1,j1) = ip
            else
               if (abs(map%imin_b(i1,j1) - im) > nibx/2) then
                  map%imin_b(i1,j1) = max(map%imin_b(i1,j1),im)
               else
                  map%imin_b(i1,j1) = min(map%imin_b(i1,j1),im)
               endif
               if (abs(map%imax_b(i1,j1) - ip) > nibx/2) then
                  map%imax_b(i1,j1) = min(map%imax_b(i1,j1),ip)
               else
                  map%imax_b(i1,j1) = max(map%imax_b(i1,j1),ip)
               endif
            endif
         enddo
         enddo

      endif

   enddo
   enddo

   !----------------------------------------------------------------------------
   write(*,F00) "fill in missing boxes and"
   write(*,F00) "make sure there are at least vcells_req valid cells per bin"
   !----------------------------------------------------------------------------

   allocate(timin(nibx,njbx),timax(nibx,njbx),tjmin(nibx,njbx),tjmax(nibx,njbx))

   cnt = 1
   cnt2 = 0
   do while (cnt > 0 .and. cnt2 < nibx)
      cnt2 = cnt2 + 1
      cnt = 0
      timin = map%imin_b
      timax = map%imax_b
      tjmin = map%jmin_b
      tjmax = map%jmax_b
      do j = 1,njbx
      do i = 1,nibx
         vcells = 0
         if (min(timin(i,j),timax(i,j),tjmin(i,j),tjmax(i,j)) > 0) then
            i1 = timin(i,j)
            i2 = timax(i,j)
            j1 = tjmin(i,j)
            j2 = tjmax(i,j)
            if (i2 > i1) then
               do j3 = j1,j2
               do i3 = i1,i2
                  n = (j3-1) * map%ni_b + i3
                  if (map%mask_b(n) > 0) vcells = vcells + 1
               enddo
               enddo
            else
               do j3 = j1,j2
               do i3 = i1,map%ni_b
                  n = (j3-1) * map%ni_b + i3
                  if (map%mask_b(n) > 0) vcells = vcells + 1
               enddo
               do i3 = 1,i2
                  n = (j3-1) * map%ni_b + i3
                  if (map%mask_b(n) > 0) vcells = vcells + 1
               enddo
               enddo
            endif
         endif
         if (vcells < vcells_req) then
            cnt = cnt + 1
            do j1 = -1,1
            do i1 = -1,1
               j2 = max(min(j+j1,njbx),1)
               i2 = mod(i+i1-1+nibx,nibx) + 1
               if (min(map%imin_b(i2,j2),map%imax_b(i2,j2),map%jmin_b(i2,j2),map%jmax_b(i2,j2)) > 0) then
                  if (min(timin(i,j),timax(i,j),tjmin(i,j),tjmax(i,j)) < 0) then
                     tjmin(i,j) = map%jmin_b(i2,j2)
                     tjmax(i,j) = map%jmax_b(i2,j2)
                     timin(i,j) = map%imin_b(i2,j2)
                     timax(i,j) = map%imax_b(i2,j2)
                  else
                     tjmin(i,j) = min(tjmin(i,j),map%jmin_b(i2,j2))
                     tjmax(i,j) = max(tjmax(i,j),map%jmax_b(i2,j2))
                     if (abs(timin(i,j)-map%imin_b(i2,j2)) > map%ni_b/2) then
                        timin(i,j) = max(timin(i,j),map%imin_b(i2,j2))
                     else
                        timin(i,j) = min(timin(i,j),map%imin_b(i2,j2))
                     endif
                     if (abs(timax(i,j)-map%imax_b(i2,j2)) > map%ni_b/2) then
                        timax(i,j) = min(timax(i,j),map%imax_b(i2,j2))
                     else
                        timax(i,j) = max(timax(i,j),map%imax_b(i2,j2))
                     endif
                  endif
               endif
            enddo
            enddo
         else

         endif

      enddo
      enddo
      map%imin_b = timin
      map%imax_b = timax
      map%jmin_b = tjmin
      map%jmax_b = tjmax
      write(6,*) subname,' fill pass number',cnt2,' count = ',cnt
   enddo

   deallocate(timin,timax,tjmin,tjmax)

!   do j = 1,njbx
!   do i = 1,nibx
!      write(6,*) 'i,j ',i,j,map%imin_b(i,j),map%imax_b(i,j),map%jmin_b(i,j),map%jmax_b(i,j)
!   enddo
!   enddo

   !----------------------------------------------------------------------------
   ! allocate matrix s,row,col, initial guess of size, will be updated
   !----------------------------------------------------------------------------

   map%n_s = map%n_a
   allocate(map%s  (map%n_s))
   allocate(map%row(map%n_s))
   allocate(map%col(map%n_s))

   do n = 1,map%n_s
      map%row(n) = n
      map%col(n) = -99
      map%S(n) = 1.0_r8
   enddo

END SUBROUTINE map_gridRead

!===============================================================================

SUBROUTINE map_gennn0(map)

   !--- modules ---

   implicit none

   !--- includes ---

   !--- arguments ---
   type(sMatrix), intent(inout) :: map       ! sMatrix info to be read in

   !--- local ---
   integer         :: i,j,n,m   ! generic indicies
   integer         :: na,nb,ns  ! index into a,b,s vectors (src,dest,matrix)
   real(r8)        :: dist      ! distance
#ifdef  _OPENMP
   integer, allocatable, dimension(:) :: nb_nn_ar ! index of nn in b (dest) vector
   real(r8), allocatable, dimension(:)    :: dmin_ar  ! minimum distance
   integer :: omp_get_thread_num, omp_get_max_threads
#endif
   integer         :: nb_nn     ! index of nn in b (dest) vector
   real(r8)            :: dmin      ! minimum distance

   integer         :: t0        ! share timer id
   character( 8)   :: cdate     ! wall clock date
   character(10)   :: ctime     ! wall clock time
   character(80)   :: str       ! 

   !--- formats ---
   character(*),parameter :: subName = "(map_gennn0) "
   character(*),parameter ::   F00 = "('(map_gennn0) ',3a)"
   character(*),parameter ::   F02 = "('(map_gennn0) ',a11,a3,60(a1))"

!-------------------------------------------------------------------------------
! NOTES:
! - create a nearest neighbor map using a "brute force" search technique
! - for every unmasked src cell there is exactly one non-zero S matrix value
!-------------------------------------------------------------------------------


   !--- allocate S: number of elements is number of unmasked src cells ---------
   map%n_s = 0
   do n = 1,map%n_a
      if (map%mask_a(n) .ne. 0) map%n_s = map%n_s + 1
   enddo
   write(6,*) ' '
   write(6,*) subname,map%n_s,' unmasked src cells out of ',map%n_a,' cells total'

   deallocate(map%s         , map%row         , map%col         )
   allocate  (map%s(map%n_s), map%row(map%n_s), map%col(map%n_s))


   !--- find nearest neighbor for each src cell --------------------------------
   call shr_timer_get  (t0,subName//"construct nearest neighbor map")
   call shr_timer_start(t0)
   ns = 0

#ifdef _OPENMP
   allocate(nb_nn_ar(omp_get_max_threads()))
   allocate(dmin_ar(omp_get_max_threads()))
#endif

   do na = 1,map%n_a
      if (mod(ns,map%n_s/50) == 0) then
         call date_and_time(cdate,ctime) ! f90 intrinsic
         str =   cdate(1:4)//'-'//cdate(5:6)//'-'//cdate(7:8)//' ' &
         &     //ctime(1:2)//':'//ctime(3:4)//':'//ctime(5:6)
         write(6,*) subname,trim(str)," done ",ns," of ",map%n_s," cells ",100*ns/map%n_s,"%"
      end if

      !--- brute force search for nearest neighbor cell ---
      dmin = 1.0e36
      nb_nn = 0
#ifdef _OPENMP
      dmin_ar(:)  = dmin
      nb_nn_ar(:) = nb_nn
#endif
      if (map%mask_a(na) .ne. 0) then ! non-zero <=> an active cell
         ns = ns + 1
         !$OMP PARALLEL DO DEFAULT(SHARED) &
         !$OMP PRIVATE(nb, dist)
         do nb = 1,map%n_b
            if (map%mask_b(nb) .ne. 0) then
               dist = map_distance(map%xc_a(na),map%yc_a(na),map%xc_b(nb),map%yc_b(nb))
#ifdef _OPENMP
               if (dist < dmin_ar(omp_get_thread_num())) then
                  dmin_ar(omp_get_thread_num())  = dist
                  nb_nn_ar(omp_get_thread_num()) = nb
               endif
#else
               if (dist < dmin) then
                  dmin = dist
                  nb_nn = nb
               endif
#endif
            endif
         enddo ! nb
         !$OMP END PARALLEL DO

#ifdef _OPENMP
         dmin = minval(dmin_ar)
         do nb=omp_get_max_threads(),1,-1
            if (dmin_ar(nb).eq.dmin) then
              nb_nn = nb_nn_ar(nb)
            end if
         end do
#endif

         if (nb_nn > 0) then
            map%col(ns) = na
            map%row(ns) = nb_nn
            map%S  (ns) = map%area_a(na)/map%area_b(nb_nn)
         else
            write(6,*) subname,'ERROR: found no nearest neighbor for src cell ',na
            call shr_sys_abort()
         endif
      endif
   end do ! na

   call shr_timer_stop (t0)
   call shr_timer_print(t0)

   if (ns /= map%n_s) then
      write(6,*) subname,' ERROR in ns ',ns,map%n_s
      call shr_sys_abort()
   endif

end subroutine map_gennn0

!===============================================================================
!===============================================================================

SUBROUTINE map_gennn(map)

   !--- modules ---

   implicit none

   !--- includes ---

   !--- arguments ---
   type(sMatrix), intent(inout) :: map       ! sMatrix info to be read in

   !--- local ---
   integer         :: i,j,n,m         ! generic indicies
   integer         :: istart,iend,jstart,jend
   integer         :: ns,cnt,ifound
   real(r8)        :: dist,dmin
   real(r8)        :: sum1,sum2

   !--- formats ---
   character(*),parameter :: subName = "(map_gennn) "
   character(*),parameter ::   F00 = "('(map_gennn) ',3a)"
   character(*),parameter ::   F02 = "('(map_gennn) ',a11,a3,60(a1))"

!-------------------------------------------------------------------------------
! NOTES:
! - needs debugging: may fail or generate nn maps with errors -kauff, 2009 Mar
!-------------------------------------------------------------------------------

   write(6,*) subname,'WARNING: this routine needs debugging.'
   call shr_sys_abort()

   deallocate(map%s  )
   deallocate(map%row)
   deallocate(map%col)

   map%n_s = 0
   do n = 1,map%n_a
      if (map%mask_a(n) > 0) map%n_s = map%n_s + 1
   enddo

   write(6,*) ' '
   write(6,*) subname,' found ',map%n_s,' points out of ',map%n_a,' points total'

   allocate(map%s  (map%n_s))
   allocate(map%row(map%n_s))
   allocate(map%col(map%n_s))

   cnt = 0
   sum1 = 0.0_r8
   sum2 = 0.0_r8
   do n = 1,map%n_a
      dmin = 1.0e36
      ifound = 0
      if (map%mask_a(n) > 0) then
         cnt = cnt + 1
         if (mod(cnt,map%n_s/50) == 1) write(6,*) subname,' doing nn ',cnt,' of ',map%n_s
         ! --- find bounding box, use +- 2 lons and lats
         call map_boundpta(map,n,istart,iend,jstart,jend)
         if (iend > istart) then
            do j = jstart,jend
            do i = istart,iend
               sum1 = sum1 + 1.0_r8
               ns = (j-1)*map%ni_b + i
               if (map%mask_b(ns) > 0) then
                  sum2 = sum2 + 1.0_r8
                  dist = map_distance(map%xc_a(n),map%yc_a(n),map%xc_b(ns),map%yc_b(ns))
                  if (dist < dmin) then
                     dmin = dist
                     ifound = ns
                  endif
               endif
            enddo
            enddo
         else
            do j = jstart,jend
            do i = 1,iend
               sum1 = sum1 + 1.0_r8
               ns = (j-1)*map%ni_b + i
               if (map%mask_b(ns) > 0) then
                  sum2 = sum2 + 1.0_r8
                  dist = map_distance(map%xc_a(n),map%yc_a(n),map%xc_b(ns),map%yc_b(ns))
                  if (dist < dmin) then
                     dmin = dist
                     ifound = ns
                  endif
               endif
            enddo
            do i = istart,map%ni_b
               sum1 = sum1 + 1.0_r8
               ns = (j-1)*map%ni_b + i
               if (map%mask_b(ns) > 0) then
                  sum2 = sum2 + 1.0_r8
                  dist = map_distance(map%xc_a(n),map%yc_a(n),map%xc_b(ns),map%yc_b(ns))
                  if (dist < dmin) then
                     dmin = dist
                     ifound = ns
                  endif
               endif
            enddo
            enddo
         endif
         if (ifound > 0) then
            map%col(cnt) = n
            map%row(cnt) = ifound
!            map%S(cnt) = map%area_b(ifound)/map%area_a(n)
            map%S(cnt) = map%area_a(n)/map%area_b(ifound)
         else
            write(6,*) subname,' ERROR in min dist search ',ifound,n,cnt
            call shr_sys_abort()
         endif
      endif
   enddo

   write(6,*) ' '
   write(6,*) subname,' avg ndst srch cells per src cell  = ',sum1/(map%n_s)
   write(6,*) subname,' percent of nonmask dst cells = ',sum2/(sum1)
   write(6,*) ' '

   if (cnt /= map%n_s) then
      write(6,*) subname,' ERROR in cnt ',cnt,map%n_s
      call shr_sys_abort()
   endif

end subroutine map_gennn

!===============================================================================

subroutine map_boundpta(map,na,istart,iend,jstart,jend)
   implicit none

   type(sMatrix),intent(in)  :: map       ! sMatrix info to write out
   integer      ,intent(in)  :: na
   integer      ,intent(out) :: istart,iend,jstart,jend

   integer :: n,im,ip,jm,jp
   integer :: nbound,ilon,ilat
   integer :: i,j,k
   integer :: nmin,dmax,num
   
   character(len=*),parameter :: subname = '(map_boundpta) '

   call map_bxindex(map%xc_a(na),map%yc_a(na),ilon,ilat)

   istart = map%imin_b(ilon,ilat)
   iend   = map%imax_b(ilon,ilat)
   jstart = map%jmin_b(ilon,ilat)
   jend   = map%jmax_b(ilon,ilat)

   if (istart < 1 .or. iend < 1 .or. istart > map%ni_b .or. iend > map%ni_b .or. &
       jstart < 1 .or. jend < 1 .or. jstart > map%nj_b .or. jend > map%nj_b) then
      write(6,*) subname,' ERROR in indices ',na,ilon,ilat,istart,iend,jstart,jend
      call shr_sys_abort()
   endif

end subroutine map_boundpta

!===============================================================================

SUBROUTINE map_print(map)

   !--- modules ---

   implicit none

   !--- includes ---

   !--- arguments ---
   type(sMatrix), intent(in)      :: map       ! sMatrix info to write out

   !--- local ---
   integer               :: i,j,n,m   ! generic indicies

   !--- formats ---
   character(len=*),parameter :: subname = '(map_print) '
   character(len=*),parameter :: F00 = "('(map_print) ',3a)"

!-------------------------------------------------------------------------------
! PURPOSE:
! o writes map information into CSM format sparse matrix data file (netCDF file)
! 
! NOTE:
!-------------------------------------------------------------------------------

   !-----------------------------------------------------------------
   ! global attributes
   !-----------------------------------------------------------------
   write(6,*) subname,' title    = ',trim(map%title)
   write(6,*) subname,' normal   = ',trim(map%normal)
   write(6,*) subname,' method   = ',trim(map%method)
   write(6,*) subname,' history  = ',trim(map%history)
   write(6,*) subname,' conventn = ',trim(map%convention)
   write(6,*) subname,' domain_a = ',trim(map%domain_a)
   write(6,*) subname,' domain_b = ',trim(map%domain_b)
   write(6,*) ' '
   write(6,*) subname,' n_a    = ',map%n_a
   write(6,*) subname,' ni_a   = ',map%ni_a
   write(6,*) subname,' nj_a   = ',map%nj_a
   write(6,*) subname,' dims_a = ',map%dims_a
   write(6,*) subname,' xc_a   = ',minval(map%xc_a),maxval(map%xc_a)
   write(6,*) subname,' yc_a   = ',minval(map%yc_a),maxval(map%yc_a)
   write(6,*) subname,' xv_a   = ',minval(map%xv_a),maxval(map%xv_a)
   write(6,*) subname,' yv_a   = ',minval(map%yv_a),maxval(map%yv_a)
   write(6,*) subname,' mask_a = ',minval(map%mask_a),maxval(map%mask_a),sum(map%mask_a)
   write(6,*) subname,' frac_a = ',minval(map%frac_a),maxval(map%frac_a),sum(map%frac_a)
   write(6,*) subname,' area_a = ',minval(map%area_a),maxval(map%area_a),sum(map%area_a)
   write(6,*) ' '
   write(6,*) subname,' n_b    = ',map%n_b
   write(6,*) subname,' ni_b   = ',map%ni_b
   write(6,*) subname,' nj_b   = ',map%nj_b
   write(6,*) subname,' dims_b = ',map%dims_b
   write(6,*) subname,' xc_b   = ',minval(map%xc_b),maxval(map%xc_b)
   write(6,*) subname,' yc_b   = ',minval(map%yc_b),maxval(map%yc_b)
   write(6,*) subname,' xv_b   = ',minval(map%xv_b),maxval(map%xv_b)
   write(6,*) subname,' yv_b   = ',minval(map%yv_b),maxval(map%yv_b)
   write(6,*) subname,' mask_b = ',minval(map%mask_b),maxval(map%mask_b),sum(map%mask_b)
   write(6,*) subname,' frac_b = ',minval(map%frac_b),maxval(map%frac_b),sum(map%frac_b)
   write(6,*) subname,' area_b = ',minval(map%area_b),maxval(map%area_b),sum(map%area_b)
   write(6,*) ' '
   write(6,*) subname,' n_s    = ',map%n_s
   write(6,*) subname,' size S = ',size(map%S),size(map%row),size(map%col)
   write(6,*) subname,' S      = ',minval(map%S),maxval(map%S)
   write(6,*) subname,' ROW    = ',minval(map%ROW),maxval(map%ROW)
   write(6,*) subname,' COL    = ',minval(map%COL),maxval(map%COL)
   write(6,*) ' '

   write(6,*) ' '
   write(6,*) subname,' Sample bin output, box number nibx,njbx = ',nibx,njbx
   do j = 1,njbx,njbx/5
   do i = 1,nibx,nibx/5
      write(6,'(a,a,2i7,a,2i7,a,2i7)') subname,' imin,imax,jmin,jmax = ',i,j,':',map%imin_b(i,j),map%imax_b(i,j),':',map%jmin_b(i,j),map%jmax_b(i,j)
   enddo
   enddo
   write(6,*) ' '

END SUBROUTINE map_print

!===============================================================================
SUBROUTINE map_write(map, filename)

   !--- modules ---

   implicit none

   !--- includes ---

   !--- arguments ---
   type(sMatrix), intent(in)      :: map       ! sMatrix info to write out
   character(*) , intent(in)      :: filename  ! name of data file

   !--- local ---
   integer               :: i,n,m   ! generic indicies
   character(len= 8)     :: cdate   ! wall clock date
   character(len=10)     :: ctime   ! wall clock time


 ! character,allocatable :: str(:)  ! variable length char string
   character(strLen)     :: str     ! variable length char string
   character(strLen)     :: attstr  ! netCDF attribute name string
   integer               :: rcode   ! netCDF routine return code
   integer               :: fid     ! netCDF file      ID
   integer               :: vid     ! netCDF variable  ID
   integer               :: did     ! netCDF dimension ID
   integer               :: vdid(2) ! vector of nc dimension ID

   !--- formats ---
   character(len=*),parameter :: F00 = "('(map_write) ',3a)"

!-------------------------------------------------------------------------------
! PURPOSE:
! o writes map information into CSM format sparse matrix data file (netCDF file)
! 
! NOTE:
!-------------------------------------------------------------------------------

   !-----------------------------------------------------------------
   ! create a new nc file
   !-----------------------------------------------------------------
!  rcode = nf_create(trim(filename),NF_CLOBBER,fid)
   rcode = nf_create(trim(filename),IOR(NF_CLOBBER,NF_64BIT_OFFSET),fid)
   if (rcode.ne.NF_NOERR) write(*,F00) nf_strerror(rcode)

   !-----------------------------------------------------------------
   ! global attributes
   !-----------------------------------------------------------------
    str  = map%title
   rcode = nf_put_att_text(fid,NF_GLOBAL,'title'      ,len_trim(str),str)
    str  = map%normal
   rcode = nf_put_att_text(fid,NF_GLOBAL,'normalization',len_trim(str),str)
    str  = map%method
   rcode = nf_put_att_text(fid,NF_GLOBAL,'map_method' ,len_trim(str),str)
    str  = "$SVN"
   rcode = nf_put_att_text(fid,NF_GLOBAL,'SVN URL'     ,len_trim(str),str)
    str  = "$Id"
   rcode = nf_put_att_text(fid,NF_GLOBAL,'SVN Id'      ,len_trim(str),str)
    str  = map%history
   call date_and_time(cdate,ctime) ! f90 intrinsic
    str = 'File created: '//cdate(1:4)//'-'//cdate(5:6)//'-'//cdate(7:8) &
   &                 //' '//ctime(1:2)//':'//ctime(3:4)//':'//ctime(5:6)
   rcode = nf_put_att_text(fid,NF_GLOBAL,'history'    ,len_trim(str),str)
    str  = map%convention
   rcode = nf_put_att_text(fid,NF_GLOBAL,'conventions',len_trim(str),str)
    str  = map%domain_a
   rcode = nf_put_att_text(fid,NF_GLOBAL,'domain_a'   ,len_trim(str),str)
    str  = map%domain_b
   rcode = nf_put_att_text(fid,NF_GLOBAL,'domain_b'   ,len_trim(str),str)

   !-----------------------------------------------------------------
   ! dimension data
   !-----------------------------------------------------------------
   rcode = nf_def_dim(fid, 'n_a' , map%n_a , did) ! # of points total
   rcode = nf_def_dim(fid, 'ni_a', map%ni_a, did) ! # of points wrt i
   rcode = nf_def_dim(fid, 'nj_a', map%nj_a, did) ! # of points wrt j
   rcode = nf_def_dim(fid, 'nv_a',  4      , did) ! # of verticies per cell
   rcode = nf_def_dim(fid, 'src_grid_rank', 2, did) ! # of verticies per cell

   rcode = nf_def_dim(fid, 'n_b' , map%n_b , did) ! # of points total
   rcode = nf_def_dim(fid, 'ni_b', map%ni_b, did) ! # of points wrt i
   rcode = nf_def_dim(fid, 'nj_b', map%nj_b, did) ! # of points wrt j
   rcode = nf_def_dim(fid, 'nv_b',  4      , did) ! # of verticies per cell
   rcode = nf_def_dim(fid, 'dst_grid_rank', 2, did) ! # of verticies per cell

   rcode = nf_def_dim(fid, 'n_s' , map%n_s , did) ! size of sparse matrix

   !-----------------------------------------------------------------
   ! define data -- coordinates, input grid
   !-----------------------------------------------------------------

   rcode = nf_inq_dimid(fid,'n_a' , did   )

   rcode = nf_def_var  (fid,'xc_a',NF_DOUBLE,1,did,vid)
   str   = 'longitude of grid cell center (input)'
   rcode = nf_put_att_text(fid,vid,"long_name",len_trim(str),str)
   str   = 'degrees east'
   rcode = nf_put_att_text(fid,vid,"units"    ,len_trim(str),str)

   rcode = nf_def_var  (fid,'yc_a',NF_DOUBLE,1,did,vid)
   str   = 'latitude of grid cell center (input)'
   rcode = nf_put_att_text(fid,vid,"long_name",len_trim(str),str)
   str   = 'degrees north'
   rcode = nf_put_att_text(fid,vid,"units"    ,len_trim(str),str)

   rcode = nf_inq_dimid(fid,'nv_a'     ,vdid(1))
   rcode = nf_inq_dimid(fid,'n_a'      ,vdid(2))

   rcode = nf_def_var  (fid,'xv_a',NF_DOUBLE,2,vdid,vid)
   str   = 'longitude of grid cell verticies (input)'
   rcode = nf_put_att_text(fid,vid,"long_name",len_trim(str),str)
   str   = 'degrees east'
   rcode = nf_put_att_text(fid,vid,"units"    ,len_trim(str),str)

   rcode = nf_def_var  (fid,'yv_a',NF_DOUBLE,2,vdid,vid)
   str   = 'latitude of grid cell verticies (input)'
   rcode = nf_put_att_text(fid,vid,"long_name",len_trim(str),str)
   str   = 'degrees north'
   rcode = nf_put_att_text(fid,vid,"units"    ,len_trim(str),str)


   rcode = nf_def_var  (fid,'mask_a',NF_INT ,1,did,vid)
   str   = 'domain mask (input)'
   rcode = nf_put_att_text(fid,vid,"long_name",len_trim(str),str)

   rcode = nf_def_var  (fid,'area_a',NF_DOUBLE ,1,did,vid)
   str   = 'area of cell (input)'
   rcode = nf_put_att_text(fid,vid,"long_name",len_trim(str),str)

   rcode = nf_def_var  (fid,'frac_a',NF_DOUBLE ,1,did,vid)
   str   = 'fraction of domain intersection (input)'
   rcode = nf_put_att_text(fid,vid,"long_name",len_trim(str),str)

   rcode = nf_inq_dimid(fid,'src_grid_rank' , did   )
   rcode = nf_def_var  (fid,'src_grid_dims',NF_INT ,1,did,vid)

   !-----------------------------------------------------------------
   ! define data -- coordinates, output grid
   !-----------------------------------------------------------------

   rcode = nf_inq_dimid(fid,'n_b'      ,did)
   if (rcode.ne.NF_NOERR) write(*,F00)  'n_b 1',nf_strerror(rcode)

   rcode = nf_def_var  (fid,'xc_b',NF_DOUBLE,1,did,vid)
   str   = 'longitude of grid cell center (output)'
   rcode = nf_put_att_text(fid,vid,"long_name",len_trim(str),str)
   str   = 'degrees east'
   rcode = nf_put_att_text(fid,vid,"units"    ,len_trim(str),str)

   rcode = nf_def_var  (fid,'yc_b',NF_DOUBLE,1,did,vid)
   str   = 'latitude of grid cell center (output)'
   rcode = nf_put_att_text(fid,vid,"long_name",len_trim(str),str)
   str   = 'degrees north'
   rcode = nf_put_att_text(fid,vid,"units"    ,len_trim(str),str)

   rcode = nf_inq_dimid(fid,'nv_b'     ,vdid(1))
   rcode = nf_inq_dimid(fid,'n_b'      ,vdid(2))

   rcode = nf_def_var  (fid,'xv_b',NF_DOUBLE,2,vdid,vid)
   str   = 'longitude of grid cell verticies (output)'
   rcode = nf_put_att_text(fid,vid,"long_name",len_trim(str),str)
   str   = 'degrees east'
   rcode = nf_put_att_text(fid,vid,"units"    ,len_trim(str),str)

   rcode = nf_def_var  (fid,'yv_b',NF_DOUBLE,2,vdid,vid)
   str   = 'latitude of grid cell verticies (output)'
   rcode = nf_put_att_text(fid,vid,"long_name",len_trim(str),str)
   str   = 'degrees north'
   rcode = nf_put_att_text(fid,vid,"units"    ,len_trim(str),str)


   rcode = nf_def_var  (fid,'mask_b',NF_INT ,1,did,vid)
   str   = 'domain mask (output)'
   rcode = nf_put_att_text(fid,vid,"long_name",len_trim(str),str)

   rcode = nf_def_var  (fid,'area_b',NF_DOUBLE ,1,did,vid)
   str   = 'area of cell (output)'
   rcode = nf_put_att_text(fid,vid,"long_name",len_trim(str),str)

   rcode = nf_def_var  (fid,'frac_b',NF_DOUBLE ,1,did,vid)
   str   = 'fraction of domain intersection (output)'
   rcode = nf_put_att_text(fid,vid,"long_name",len_trim(str),str)

   rcode = nf_inq_dimid(fid,'dst_grid_rank' , did   )
   rcode = nf_def_var  (fid,'dst_grid_dims',NF_INT ,1,did,vid)

   !-----------------------------------------------------------------
   ! define data -- matrix elements
   !-----------------------------------------------------------------

   !code = nf_inq_dimid(fid,'n_wgts',vdid(1)) ! no gradients!
   rcode = nf_inq_dimid(fid,'n_s'   , did   )
   rcode = nf_def_var  (fid,'S',NF_DOUBLE,1,did,vid)
   str   = 'sparse matrix for mapping S:a->b'
   rcode = nf_put_att_text(fid,vid,"long_name",len_trim(str),str)

   rcode = nf_inq_dimid(fid,'n_s'   ,did)

   rcode = nf_def_var  (fid,'col',NF_INT ,1,did,vid)
   str   = 'column corresponding to matrix elements'
   rcode = nf_put_att_text(fid,vid,"long_name",len_trim(str),str)

   rcode = nf_def_var  (fid,'row',NF_INT ,1,did,vid)
   str   = 'row corresponding to matrix elements'
   rcode = nf_put_att_text(fid,vid,"long_name",len_trim(str),str)

   !-----------------------------------------------------------------
   ! put data
   !-----------------------------------------------------------------
   rcode = nf_enddef(fid)

   rcode = nf_inq_varid     (fid,'src_grid_dims', vid)
   rcode = nf_put_var_int   (fid,    vid , map%dims_a)
   rcode = nf_inq_varid     (fid,'xc_a', vid)
   rcode = nf_put_var_double(fid,  vid ,map%xc_a)
   rcode = nf_inq_varid     (fid,'yc_a', vid)
   rcode = nf_put_var_double(fid,  vid ,map%yc_a)
   rcode = nf_inq_varid     (fid,'xv_a', vid)
   rcode = nf_put_var_double(fid,  vid ,map%xv_a)
   rcode = nf_inq_varid     (fid,'yv_a', vid)
   rcode = nf_put_var_double(fid,  vid ,map%yv_a)
   rcode = nf_inq_varid     (fid,'mask_a',vid)
   rcode = nf_put_var_int   (fid,    vid , map%mask_a)
   rcode = nf_inq_varid     (fid,'area_a', vid)
   rcode = nf_put_var_double(fid,  vid ,map%area_a)
   rcode = nf_inq_varid     (fid,'frac_a', vid)
   rcode = nf_put_var_double(fid,  vid ,map%frac_a)

   rcode = nf_inq_varid     (fid,'dst_grid_dims', vid)
   rcode = nf_put_var_int   (fid,    vid , map%dims_b)
   rcode = nf_inq_varid     (fid,'xc_b', vid)
   rcode = nf_put_var_double(fid,  vid ,map%xc_b)
   rcode = nf_inq_varid     (fid,'yc_b', vid)
   rcode = nf_put_var_double(fid,  vid ,map%yc_b)
   rcode = nf_inq_varid     (fid,'xv_b', vid)
   rcode = nf_put_var_double(fid,  vid ,map%xv_b)
   rcode = nf_inq_varid     (fid,'yv_b', vid)
   rcode = nf_put_var_double(fid,  vid ,map%yv_b)
   rcode = nf_inq_varid     (fid,'mask_b',vid)
   rcode = nf_put_var_int   (fid,    vid , map%mask_b)
   rcode = nf_inq_varid     (fid,'area_b', vid)
   rcode = nf_put_var_double(fid,  vid ,map%area_b)
   rcode = nf_inq_varid     (fid,'frac_b', vid)
   rcode = nf_put_var_double(fid,  vid ,map%frac_b)

   rcode = nf_inq_varid     (fid,  'S',vid)
   rcode = nf_put_var_double(fid, vid, map%s) ! no gradients!
   rcode = nf_inq_varid     (fid,'row',vid)
   rcode = nf_put_var_int   (fid, vid, map%row)
   rcode = nf_inq_varid     (fid,'col',vid)
   rcode = nf_put_var_int   (fid, vid, map%col)

   rcode = nf_close(fid)

   if (rcode.ne.NF_NOERR) write(*,F00) nf_strerror(rcode)


END SUBROUTINE map_write

!===============================================================================

SUBROUTINE map_Sn1(map)

   implicit none

   !--- arguments ---
   type(sMatrix),intent(inout) :: map

   !--- local ---
   integer :: i,n

!-------------------------------------------------------------------------------
! PURPOSE:
! X compute sn1(i) = least n st row(n)=i 
!                  = least link number specifying a mapping
!                    to the ith destination point (ith row number).
! o compute sn1(i) = # of links for this row
! o compute sn2(i) = total # of links for row < i
!-------------------------------------------------------------------------------

!  i=0 ! last row for which sn1(i) was calculated
!  do n=1,map%n_s
!    if (map%row(n) /= i) then
!      map%sn1(i+1:map%row(n)) = n
!      i=map%row(n)
!    endif
!  enddo
!  map%sn1(i+1:map%n_b+1)=map%n_s+1

   map%sn1 = 0
   map%sn2 = 0
   do n=1,map%n_s
     map%sn1(map%row(n)) = map%sn1(map%row(n))+1
   enddo
   do n=2,map%n_b
     map%sn2(n) = map%sn2(n-1) + map%sn1(n-1)
   enddo

END SUBROUTINE map_Sn1

!===============================================================================

SUBROUTINE map_dup(map_in,map_out)

   implicit none

   !--- arguments ---
   type(sMatrix),intent( in)   :: map_in
   type(sMatrix),intent(inout) :: map_out

   !--- local ---
   integer :: rcode ! return code

!-------------------------------------------------------------------------------
! PURPOSE:
! o create a duplicate copy a sparse_matrix (cannot deallocate existing data)
!
! NOTES
! o the return code allows the deallocate to fail without stopping program 
!   exectution.  This would happen, eg, if the pointer had never been 
!   allocated to begin with (dealloc doesn't like undefined pointers). 
!   F90 has no way (?) of testing an undefined pointer (eg. never allocated).
!   See Fortran 90/95 Explained, Melcaf, 1998, p 124, 175
!-------------------------------------------------------------------------------

   !------------------------------------------------
   ! de-allocate space
   !------------------------------------------------
   deallocate(map_out%  xc_a,STAT=rcode) 
   deallocate(map_out%  yc_a,STAT=rcode) 
   deallocate(map_out%  xv_a,STAT=rcode) 
   deallocate(map_out%  yv_a,STAT=rcode) 
   deallocate(map_out%mask_a,STAT=rcode) 
   deallocate(map_out%area_a,STAT=rcode) 

   deallocate(map_out%  xc_b,STAT=rcode) 
   deallocate(map_out%  yc_b,STAT=rcode) 
   deallocate(map_out%  xv_b,STAT=rcode) 
   deallocate(map_out%  yv_b,STAT=rcode) 
   deallocate(map_out%mask_b,STAT=rcode) 
   deallocate(map_out%area_b,STAT=rcode) 

   deallocate(map_out%frac_a,STAT=rcode) 
   deallocate(map_out%frac_b,STAT=rcode) 

   deallocate(map_out%s     ,STAT=rcode) 
   deallocate(map_out%row   ,STAT=rcode) 
   deallocate(map_out%col   ,STAT=rcode) 
   deallocate(map_out%sn1   ,STAT=rcode) 
   deallocate(map_out%sn2   ,STAT=rcode) 

   !------------------------------------------------
   ! allocate space
   !------------------------------------------------
   allocate(map_out%  xc_a(  map_in%n_a) )
   allocate(map_out%  yc_a(  map_in%n_a) )
   allocate(map_out%  xv_a(4,map_in%n_a) )
   allocate(map_out%  yv_a(4,map_in%n_a) )
   allocate(map_out%mask_a(  map_in%n_a) )
   allocate(map_out%area_a(  map_in%n_a) )

   allocate(map_out%  xc_b(  map_in%n_b) )
   allocate(map_out%  yc_b(  map_in%n_b) )
   allocate(map_out%  xv_b(4,map_in%n_b) )
   allocate(map_out%  yv_b(4,map_in%n_b) )
   allocate(map_out%mask_b(  map_in%n_b) )
   allocate(map_out%area_b(  map_in%n_b) )

   allocate(map_out%frac_a(map_in%n_a) )
   allocate(map_out%frac_b(map_in%n_b) )  

   allocate(map_out%s     (map_in%n_s) )  
   allocate(map_out%row   (map_in%n_s) )  
   allocate(map_out%col   (map_in%n_s) )  
   allocate(map_out%sn1   (map_in%n_b) )  
   allocate(map_out%sn2   (map_in%n_b) )  

   !------------------------------------------------
   ! set values
   !------------------------------------------------
   map_out%   n_a = map_in%   n_a
   map_out%  ni_a = map_in%  ni_a
   map_out%  nj_a = map_in%  nj_a
   map_out%dims_a = map_in%dims_a
   map_out%  xc_a = map_in%  xc_a
   map_out%  yc_a = map_in%  yc_a
   map_out%  xv_a = map_in%  xv_a
   map_out%  yv_a = map_in%  yv_a
   map_out%mask_a = map_in%mask_a
   map_out%area_a = map_in%area_a

   map_out%   n_b = map_in%   n_b
   map_out%  ni_b = map_in%  ni_b
   map_out%  nj_b = map_in%  nj_b
   map_out%dims_b = map_in%dims_b
   map_out%  xc_b = map_in%  xc_b
   map_out%  yc_b = map_in%  yc_b
   map_out%  xv_b = map_in%  xv_b
   map_out%  yv_b = map_in%  yv_b
   map_out%mask_b = map_in%mask_b
   map_out%area_b = map_in%area_b

   map_out%frac_a = map_in%frac_a
   map_out%frac_b = map_in%frac_b

   map_out%n_s    = map_in%n_s
   map_out%s      = map_in%s
   map_out%row    = map_in%row
   map_out%col    = map_in%col
   map_out%sn1    = map_in%sn1
   map_out%sn2    = map_in%sn2

   map_out%title      = map_in%title
   map_out%normal     = map_in%normal
   map_out%method     = map_in%method
   map_out%history    = map_in%history
   map_out%convention = map_in%convention
   map_out%domain_a   = map_in%domain_a
   map_out%domain_b   = map_in%domain_b

END SUBROUTINE map_dup

!===============================================================================

SUBROUTINE map_free(map)

   !--- modules ---

   implicit none

   !--- includes ---

   !--- arguments ---
   type(sMatrix), intent(inout) :: map       ! sMatrix to be deallocated

   !--- local ---
   integer :: rcode ! return code

   !--- formats ---
   character(len=*),parameter :: F00 = "('(map_free) ',3a)"

!-------------------------------------------------------------------------------
! PURPOSE:
! o deallocate a map matrix data type
! 
! NOTES
! o the return code allows the deallocate to fail without stopping program 
!   exectution.  This would happen, eg, if the pointer had never been 
!   allocated to begin with (dealloc doesn't like undefined pointers). 
!   F90 has no way (?) of testing an undefined pointer (eg. never allocated).
!   See Fortran 90/95 Explained, Melcaf, 1998, p 124, 175
!-------------------------------------------------------------------------------

   !-----------------------------------------------
   ! unset global attributes
   !-----------------------------------------------
   map%title      = "null-free"
   map%normal     = "null-free"
   map%method     = "null-free"
   map%convention = "null-free"
   map%history    = "null-free"
   map%domain_a   = "null-free"
   map%domain_b   = "null-free"

   !-----------------------------------------------
   ! free "a" domain info 
   !-----------------------------------------------

   deallocate(map%  xc_a,STAT=rcode)
   deallocate(map%  yc_a,STAT=rcode)
   deallocate(map%  xv_a,STAT=rcode)
   deallocate(map%  yv_a,STAT=rcode)
   deallocate(map%mask_a,STAT=rcode)
   deallocate(map%area_a,STAT=rcode)
   deallocate(map%frac_a,STAT=rcode)

   !-----------------------------------------------
   ! free "b" domain info 
   !-----------------------------------------------

   deallocate(map%  xc_b,STAT=rcode)
   deallocate(map%  yc_b,STAT=rcode)
   deallocate(map%  xv_b,STAT=rcode)
   deallocate(map%  yv_b,STAT=rcode)
   deallocate(map%mask_b,STAT=rcode)
   deallocate(map%area_b,STAT=rcode)
   deallocate(map%frac_b,STAT=rcode)

   !-----------------------------------------------
   ! free matrix info 
   !-----------------------------------------------
   deallocate(map%s  ,STAT=rcode)
   deallocate(map%row,STAT=rcode)
   deallocate(map%col,STAT=rcode)
   deallocate(map%sn1,STAT=rcode)
   deallocate(map%sn2,STAT=rcode)

END SUBROUTINE map_free

!===============================================================================

subroutine map_bxindex(lon,lat,ibx,jbx)

   implicit none

   !--- arguments ---
   real(r8)   , intent(in) :: lon,lat
   integer, intent(out) :: ibx,jbx

   !--- local ---
   real(r8) :: llon,llat
   character(*),parameter :: subName = "(map_bxindex) "
   real(r8),parameter :: lonmin =  0.0_r8
   real(r8),parameter :: lonmax =  360.0_r8
   real(r8),parameter :: latmin = -90.0_r8
   real(r8),parameter :: latmax =  90.0_r8


!-------------------------------------------------------------------------------
! PURPOSE:
!   derives the bx index for a given lon, lat.
!-------------------------------------------------------------------------------

   llat = lat
   llon = mod(lon+720.0_r8,360.0_r8)

   if (llat < latmin .or. llat > latmax) then
      write(6,*) subname,' ERROR in lat ',llat
      call shr_sys_abort()
   endif

   if (llon < lonmin .or. llon > lonmax) then
      write(6,*) subname,' ERROR in lon ',llon
      call shr_sys_abort()
   endif

   if (llon == lonmax) then
      ibx = nibx
   else
      ibx = ((llon - lonmin) / (lonmax -lonmin)) * (nibx * 1.0_r8) + 1
   endif
   if (llat == latmax) then
      jbx = njbx
   else
      jbx = ((llat - latmin) / (latmax -latmin)) * (njbx * 1.0_r8) + 1
   endif

   if (ibx < 1 .or. ibx > nibx .or. jbx < 1 .or. jbx > njbx) then
      write(6,*) subname,' ERROR in ibx,jbx ',ibx,jbx,nibx,njbx
   endif
   

end subroutine map_bxindex

!===============================================================================
!===============================================================================

integer FUNCTION nSij(smat,i,j)

   implicit none

   !--- arguments ---
   type(sMatrix),intent(in ) :: smat ! sparse matrix data type
   integer      ,intent(in ) :: i,j  ! index into A(i,j)

   !--- local ---
   integer    :: n                   ! generic index 

!-------------------------------------------------------------------------------
! PURPOSE:
!   find n st (i,j) = (row(n),col(n))
!   else return 0
!-------------------------------------------------------------------------------

!  nSij=0
!  do n=smat%sn1(i),smat%sn1(i+1)-1
!    if ( smat%col(n) == j ) exit
!  end do
!  if ( smat%row(n)/=i   .or.  smat%col(n)/=j ) then
!     nSij = 0
!  else
!     nSij = n
!  end if

   nSij=0
   if (smat%sn1(i) > 0) then
      do n=smat%sn2(i)+1,smat%sn2(i)+smat%sn1(i)
        if (smat%col(n) == j) exit
      enddo
      if ( smat%row(n)/=i   .or.  smat%col(n)/=j ) then
         nSij = 0
      else
         nSij = n
      end if
   endif

END FUNCTION nSij

!===============================================================================

SUBROUTINE map_matMatMult(A,B,S)

   implicit none

   !--- arguments ---
   type(sMatrix),intent(in )   :: A  ! sparse matrix data type: B=SA
   type(sMatrix),intent(in )   :: S  ! sparse matrix data type: B=SA
   type(sMatrix),intent(inout) :: B  ! sparse matrix data type: B=SA

   !--- local ---
   integer    :: ib,jb             ! index wrt B(j,i)
   integer    :: ia,ja             ! index wrt A(j,i)
   integer    :: is,js             ! index wrt S(i,j)
   integer    :: i,j,k             ! index wrt A(j,k),S(i,j),B(j,k)
   integer    :: ns                ! index wrt S(n) 
   integer    :: m                 ! index wrt b(i,m),a(j,m)
   integer    :: n                 ! index wrt S%S(n)
   integer    :: na                ! index wrt A%S(n)
   integer    :: nb                ! index wrt B%S(n)
   integer    :: nbMax             ! max size of B%S(n)
   real(r8)   :: sum               ! register variable for b(i,m)
   integer    :: t1 = 0            ! system-clock-timer number
   integer    :: t0 = 0            ! system-clock-timer number
   integer    :: nMod              ! for periodic timer output
   integer    :: ni
   integer,pointer :: itemp(:)
   real(r8)   ,pointer :: temp(:)
   integer    :: esize

#ifdef  _OPENMP
   integer  :: omp_get_max_threads !  $OMP function call
#endif

   !--- formats ---
   character(len=*),parameter :: F00 = "('<map_matMatMult> ',a,50(f6.2 ))"
   character(len=*),parameter :: F01 = "('<map_matMatMult> ',a,50(i3,3x))"
   character(len=*),parameter :: F02 = "('<map_matMatMult> ',    f7.3 )"
   character(len=*),parameter :: F03 = "('<map_matMatMult> ',a,i5,a,i2)"
   character(len=*),parameter :: F04 = "('<map_matMatMult> ',a,2i10)"

   character(8)   :: dstr       ! F90 wall clock date string yyyymmdd
   character(10)  :: tstr       ! F90 wall clock time string hhmmss.sss
   character(*),parameter :: F10 = "('(map_matMatMult) ',a,i7,' date & time:',1x, &
                                 &   a4,2('-',a2),2x,a2,2(':',a2))"
   character(*),parameter :: subName = "(map_matMatMult) "

!-------------------------------------------------------------------------------
! PURPOSE:
!   does a sparse-matrix-multiply: B = S*A 
! ASSUMPTION:
!   o matricies are st for S(n) ~ S(i,j), 
!     * i(n) is monotonically increasing
!     * if i(n) = constant for n = n0,..n1 
!       then j(n) is monotonically increasing for n = n0,...n1
!   o B%s(:) is sufficiently large to store incoming data
!-------------------------------------------------------------------------------


   if (t0 == 0) then
      call shr_timer_get  (t0,subName//"matrix multiply - all")
      call shr_timer_get  (t1,subName//"matrix multiply - 100k rows")
      nMod = 100000 ! "100k" referred to by timer
#ifdef  _OPENMP
     n = omp_get_max_threads()
     write(6,F03) 'FYI: this routine is threaded, omp_get_max_threads() = ',n
#else
     write(6,F03) 'FYI: this routine is not threaded'
#endif
   end if

   if (B%n_a /= A%n_a) then
     write(6,*) 'ERROR dim 2 of A & B not the same'
     stop
   end if
   if (B%n_b /= S%n_b) then
     write(6,*) 'ERROR dim 1 of B & S not the same'
     stop
   end if
   if (S%n_a /= A%n_b) then
     write(6,*) 'ERROR dim 2 of S & dim 1 A not the same'
     stop
   end if

   !-----------------------------------------------------------------
   ! matrix multiply 
   !-----------------------------------------------------------------

   write(6,F04) 'matrix-matrix multiply: B = S*A'
   write(6,F04) 'do i=1 to number rows in matrix S:',S%n_b
   write(6,F04) 'do i=k to number cols in matrix A:',A%n_a
   call shr_timer_start(t0)

   !--- array size estimate ---
   esize = (S%n_s/min(S%n_b,A%n_s)) * A%n_s * 2
   esize = int(float(esize) * 1.2)  ! this was necessary for some paleo mapping file - kauff

   deallocate( B%s  )
   deallocate( B%row)
   deallocate( B%col)
   allocate(B%s  (esize))
   allocate(B%row(esize))
   allocate(B%col(esize))
   B%n_s = 0

#if (1 == 0)
   !--- assume A has 1 value per column -- implementation one (slow)

   n = 0
   do i = 1,S%n_b
      do k = 1,A%n_s
         do ni=S%sn2(i)+1,S%sn2(i)+S%sn1(i)
            if (A%row(k) == S%col(ni)) then
               n = n + 1
               if (n > esize) then
                  write(6,*) subname,' ERROR esize too small',esize,n
                  call shr_sys_abort()
               endif
               B%row(n) = S%row(ni)
               B%col(n) = A%col(k)
               B%S(n) = A%S(k) * S%S(ni)
            endif
         enddo
      enddo
      if (mod(i,S%n_b/50) == 0) write(6,*) subname,' done ',i,' of ',S%n_b,'....',100*i/S%n_b,' percent, nwgts = ',n
      if (mod(n,1000000) == 0) write(6,*) subname,' done ',i,' of ',S%n_b,'....',100*i/S%n_b,' percent, nwgts = ',n
   enddo
   B%n_s = n

#endif

#if (1 == 1)
   !--- assume A has 1 value per column -- implementation two

   n = 0
   do i = 1,S%n_b
      do ni=S%sn2(i)+1,S%sn2(i)+S%sn1(i)
         j = S%col(ni)
!         if (S%row(ni) /= i) then
!            write(6,*) subname,' ERROR in s%row(ni),i,ni = ',s%row(ni),i,ni
!            call shr_sys_abort()
!         endif
         do k = A%sn2(j)+1,A%sn2(j)+A%sn1(j)
            n = n + 1
            if (n > esize) then
               write(6,*) subname,' ERROR esize too small',esize,n
               call shr_sys_abort()
            endif
            B%row(n) = S%row(ni)
            B%col(n) = A%col(k)
            B%S(n) = S%S(ni) * A%S(k)
         enddo
      enddo
      if (mod(i,S%n_b/50) == 0) write(6,*) subname,' done ',i,' of ',S%n_b,'....',100*i/S%n_b,' percent'
   enddo
   B%n_s = n

#endif

#if (1 == 0)
   !--- original implementation for general matrix, now hangs with NN maps?

   nbMax = size(B%s)
   nb=0 ! number of non-zero elements in B
!  MP PARALLEL DO PRIVATE(i,j,k,n,na,sum)
   do i=1,S%n_b    !--- i loops over rows of S & B ---
      !------------------------------------------------------
      ! progress report
      !------------------------------------------------------
      if (mod(i-1,nMod)==0) then
         if (i > 1) then
            call shr_timer_stop (t1)
            call shr_timer_print(t1)
            call shr_sys_flush(6)
         end if
         call shr_timer_zero (t1)
         call shr_timer_start(t1)
       ! if (mod(i-1,10000)==0) then
            call date_and_time(dstr,tstr) ! F90 intrinsic routine
            write(6,F10) 'i=',i-1,dstr(1:4),dstr(5:6),dstr(7:8),tstr(1:2),tstr(3:4),tstr(5:6)
       ! end if
      end if
!$OMP PARALLEL DO PRIVATE(j,k,n,na,sum)
      do k=1,A%n_a  !--- k loops over cols of A & B ---
         !------------------------------------------------------
         ! compute B(i,k)
         !------------------------------------------------------
         sum = 0.0
         do n=S%sn2(i)+1,S%sn2(i)+S%sn1(i)  !--- j loops over S cols & A rows ---
!            write(6,*) 'tcx matmult1 ',i,n,s%sn2(i)+1,s%sn2(i)+S%sn1(i)
            j  = S%col(n)          !--- need  S(i,j)*A(j,k) for j=...
            na = nSij(A,j,k)       !--- A(na) = A(j,k)
!            write(6,*) 'tcx matmult2 ',n,j,k,na,S%S(n),A%S(na)
            if (na>0) sum = sum + S%S(n)*A%S(na)
         end do
         !------------------------------------------------------
         ! if B(i,k) > 0, add element to sparse matrix data
         !------------------------------------------------------
         if (sum /= 0.0) then 
!$OMP CRITICAL
            nb    = nb+1
            B%n_s = nb
            n     = nb
!$OMP END CRITICAL
            if (n > nbMax ) then
               write(6,*) 'ERROR: nb > nbMax (fatal)'
               call shr_sys_abort(subName//"ERROR: nb > nbMax")
            !  stop ! OMP has a problem with this
            end if
            B%S  (n) = sum
            B%row(n) = i
            B%col(n) = k
         end if
         !------------------------------------------------------
      end do
!$OMP BARRIER
   end do

#endif

   write(6,*) ' '
   write(6,*) subname,' computed ',B%n_s,' weights, size is ',size(B%s)
   write(6,*) ' '
                 
   n = B%n_s 
   allocate(itemp(n))

   itemp(1:n) = B%row(1:n)
   deallocate(B%row)
   allocate(B%row(n))
   B%row(1:n)  = itemp(1:n)

   itemp(1:n) = B%col(1:n)
   deallocate(B%col)
   allocate(B%col(n))
   B%col(1:n)  = itemp(1:n)

   deallocate(itemp)
   allocate(temp(n))

   temp(1:n) = B%S(1:n)
   deallocate(B%S)
   allocate(B%S(n))
   B%S(1:n)  = temp(1:n)

   deallocate(temp)

   call shr_timer_stop (t0)
   call shr_timer_print(t0)

END SUBROUTINE map_matMatMult

!===============================================================================

SUBROUTINE map_check(sMat)

   implicit none

   !--- arguments ---
   type(sMatrix),intent(in )   :: sMat    ! sparse matrix data type: B=SA

   !--- local ---
   integer        :: i,j,k,m,n,p          ! generic index
   integer        :: n0,n1,n2,n3,n4,n5,n6 ! generic counters
   integer        :: n7,n8,n9             ! generic counters
   integer        :: oldrow               ! generic counter
   integer        :: thisrow, thiscol     ! generic counter
   real(r8)       :: sum,maxsum,minsum    ! generic summations
   real(r8)       :: maxerr,avgerr
   real(r8)       :: f                    ! temporary float variable
   real(r8)       :: dist,maxdist         ! distance, max distance
   real(r8)       :: dph,dth              ! dphi,dtheta wrt cell distances
   real(r8)       :: x0,x1,y0,y1,dx,dy    ! 
   real(r8)       :: mn,mx                ! min & max
   integer        :: ncol,mincol,maxcol   ! # of cols per row
   real(r8),parameter :: eps = 1.0e-1         ! epsilon wrt X =? 0.0
   real(r8),allocatable :: colsum(:)          ! column sum w/o  dest area fractions
   real(r8),allocatable :: colfsum(:)         ! column sum with dest area fractions
   integer,allocatable :: colnsum(:)      ! column sum with dest area fractions

   !--- formats ---
   character(len=*),parameter :: F00 = "('(map_check) ',a,5i11)"
   !----------------------------------------------------------------------------
   character(len=*),parameter :: F01 = "('(map_check) row sum    ', &
   &  '     < 0        ~ 0      (0,1)        ~ 1        > 1     epsilon')"
   character(len=*),parameter :: F02 = "('(map_check) row sum    ', &
   &  '---------- ---------- ---------- ---------- ---------- ----------')"
   character(len=*),parameter :: F03 = "('(map_check) row sum # ',5i11,es11.3)"
   character(len=*),parameter :: F04 = "('(map_check) row sum % ',5F11.6)"
   character(len=*),parameter :: F05 = "('(map_check) row sum min/max',2es11.3)"
   character(len=*),parameter :: F06 = "('(map_check) row sum ',&
   & 'min/max cols per row ',2i6)"
   !----------------------------------------------------------------------------
   character(len=*),parameter :: f08 = "('(map_check) col(:) min/max',2i11)"
   character(len=*),parameter :: F08a = "('(map_check) frac_a(:) min/max',2f15.10)"
   character(len=*),parameter :: F09 = "('(map_check) row(:) min/max',2i11)"
   character(len=*),parameter :: F09a = "('(map_check) frac_b(:) min/max',2f15.10)"
   character(len=*),parameter :: F10 = "('(map_check) S(:)   min/max',2es11.3)"
   !----------------------------------------------------------------------------
   character(len=*),parameter :: F11 = "('(map_check) ',a,i11,2a)"
   !----------------------------------------------------------------------------
   character(len=*),parameter :: F12 ="('(map_check) ',a,f15.12,'  basis case')"
   character(len=*),parameter :: f13 = "('(map_check) ',a,f15.12,f13.6,'%')"
   !----------------------------------------------------------------------------
   character(len=*),parameter :: F20 = "('(map_check) ',5a)"
   character(len=*),parameter :: F21 = "('(map_check) ',a,i11,f13.6,'%')"
   character(len=*),parameter :: F22 = "('(map_check) ',a,2i11)"
   character(len=*),parameter :: F23 = "('(map_check) ',a,2f13.6)"

!-------------------------------------------------------------------------------
! PURPOSE:
!   perform a standard battery of tests on a matrix file
!-------------------------------------------------------------------------------

   allocate(colsum (sMat%n_a))
   allocate(colfsum(sMat%n_a))
   allocate(colnsum(sMat%n_a))

   write(6,F20) '=============================================================='
   write(6,F20) 'matrix title: ',trim(sMat%title)
   write(6,F00) 'n_s,n_a,n_b=',sMat%n_s,sMat%n_a,sMat%n_b


if (.true.) then
   write(6,F20) '--------------------------------------------------------------'
   write(6,F00) 'check for correct row/column sorting                          '
   write(6,F20) '--------------------------------------------------------------'
   i=0
   do n=2,sMat%n_s
     if (sMat%row(n-1) > sMat%row(n) ) i=i+1
   enddo
   if (i==0) then
      write(6,F00)    "OK: row(:) is monotonically increasing"
   else
      write(6,F00) "ERROR: row(:) is not monotonically increasing"
   end if

   j=0
   do i=1,sMat%n_s-1
      if (sMat%col(i).gt.sMat%col(i+1).and.sMat%row(i+1).eq.sMat%row(i)) then
         j=j+1
      endif
   enddo
   if (j==0) then
      write(6,F00)    "OK: for row=i, col(:) is monotonically increasing"
   else
      write(6,F00) "WARNING! for row=i, col(:) is not monotonically increasing"
   end if

   !----------------------------------------------------------------------------
   write(6,F00) 'Do S,row,col contain reasonable values?                       '
   !----------------------------------------------------------------------------

   write(6,F08) minval(sMat%col(1:sMat%n_s)),maxval(sMat%col(1:sMat%n_s))
   write(6,F08a) minval(sMat%frac_a(1:sMat%n_a)),maxval(sMat%frac_a(1:sMat%n_a))
   write(6,F09) minval(sMat%row(1:sMat%n_s)),maxval(sMat%row(1:sMat%n_s))
   write(6,F09a) minval(sMat%frac_b(1:sMat%n_b)),maxval(sMat%frac_b(1:sMat%n_b))
!  write(6,F10) minval(sMat%s  (1:sMat%n_s)),maxval(sMat%s  (1:sMat%n_s))
!  too big for Blackforest stack size???
   mn =  1.0e30
   mx = -1.0e30
   do n=1,sMat%n_s
     mn = min(mn,sMat%s(n))
     mx = max(mx,sMat%s(n))
   end do
   write(6,F10) mn,mx

   !----------------------------------------------------------------------------
   write(6,F00) 'Are col->row links unique? (assumes mapping is sorted)    '
   !----------------------------------------------------------------------------

   thisrow = 0
   n1 = 0
   do n=1,sMat%n_s
     if (sMat%row(n) /= thisrow) then
       thisrow = sMat%row(n)
       do m=0,sMat%sn1(thisrow)-2
         if (sMat%col(n+m) == sMat%col(n+m+1)) n1 = n1 + 1
         if (sMat%row(n+m) /= sMat%row(n+m+1)) then
!        write(6,F22) " ERROR:  ",thisrow,sMat%sn1(thisrow)
         endif
       enddo
     endif
   enddo
   if (n1 > 0) then 
     write(6,F11) "WARNING!  found nonunique links between col->row  ",n1
   else
     write(6,F11) "OK, this matrix has unique col->row links"
   endif

   write(6,F20) '--------------------------------------------------------------'
   write(6,F00) 'Definitive conservation check                                 '
   write(6,F20) '--------------------------------------------------------------'

   colsum =0.0
   colnsum=0
   do n=1,sMat%n_s
       i = sMat%row(n)
       j = sMat%col(n)
        colsum(j) =  colsum(j) + sMat%s(n)*sMat%area_b(i)
       colnsum(j) = colnsum(j) + 1
   enddo
   n0=0
   n1=0
   n2=0
   n3=0
   n4=0
   n5=0
   n6=0
   n7=0
   n8=0
   n9=0
   maxerr = 0.0
   avgerr = 0.0
   do j=1,sMat%n_a
     if ( sMat%mask_a(j) /= 0) then
       n0 = n0 + 1
       if ( colnsum(j)==0 ) n5 = n5+1
       if ( sMat%frac_a(j) == 0.0 ) n8 = n8+1
       if ( sMat%frac_a(j) > 0.0 .and. sMat%frac_a(j) < 1.0 ) n9 = n9+1
   !!! if ( abs( colsum(j)-sMat%area_a(j) )/sMat%area_a(j)  < .0001) then
   !   if ( abs( colsum(j)-sMat%area_a(j) )/sMat%area_a(j)  < .001 ) then
   !   if ( abs( colsum(j)-sMat%area_a(j) )/sMat%area_a(j)  < .01 ) then
       if ( abs( colsum(j)-sMat%area_a(j)*sMat%frac_a(j) )/ &
  &               sMat%area_a(j)*sMat%frac_a(j)  < .0000001 ) then
          n1=n1+1
       else
          n2=n2+1
          avgerr = avgerr + abs( colsum(j)-sMat%area_a(j) )/sMat%area_a(j)
          if ( abs( colsum(j)-sMat%area_a(j) )/sMat%area_a(j) > maxerr) then
            maxerr = abs( colsum(j)-sMat%area_a(j) )/sMat%area_a(j)
          end if 
       end if
     else
       if ( colnsum(j) == 0.0) then
          n3=n3+1
       else
          n4=n4+1
       end if
     end if
   enddo
   oldrow = 0
   do j=1,sMat%n_s
     if (sMat%row(j) /= oldrow) then
       if ( sMat%mask_b(sMat%row(j)) /= 0) then
         n6 = n6 + 1
       else
         n7 = n7+1
       end if
       oldrow = sMat%row(j)
     end if
   enddo
   write(6,F22) ' "correct column sum" means' 
   write(6,F22) '  sum over i of S(i,j)*area_b(i) = area_a(j)'
   write(6,F22) "  number of columns                              ",sMat%n_a
   write(6,F22) "  number of active columns                       ",n0
   write(6,F22) "  number of active columns with correct sum      ",n1
   write(6,F22) "  number of active columns with incorrect sum    ",n2
   write(6,F23) "  avg % error of incorrect sum                   ",avgerr/(max(n2,1)*1.0_r8)
   write(6,F23) "  max % error of incorrect sum                   ",maxerr
   write(6,F22) "  number of active columns with no links         ",n5
   write(6,F22) "  number of active columns with frac_a = 0.0     ",n8
   write(6,F22) "  number of active columns with 0.0<frac_a<1.0   ",n9
   write(6,F22) "  number of inactive columns                     ",sMat%n_a-n0
   write(6,F22) "  number of inactive columns with no links       ",n3
   write(6,F22) "  number of inactive columns with some links     ",n4
   write(6,F22) "  number of rows                                 ",sMat%n_b
   write(6,F22) "  number of inactive rows with some links        ",n7
end if

if (.false.) then
   !----------------------------------------------------------------------------
   write(6,F00) 'mapping distance...'
   !----------------------------------------------------------------------------
   maxdist=0.0
   sum    =0.0
   dx = 0.0 ! max dx distance
   dy = 0.0 ! max dy distance
   
   do n=1,sMat%n_s
       i = sMat%row(n)
       j = sMat%col(n)

       x0 = sMat%xc_a(j)
       x1 = sMat%xc_b(i)
       y0 = sMat%yc_a(j)
       y1 = sMat%yc_b(i)

       if (y1<95.0 .and. y0<95.0) then

         if (x1-x0>180.0) x0 = x0 + 360.0
         if (x0-x1>180.0) x1 = x1 + 360.0

         dx = max(dx,abs(x1-x0)*cos(0.5*(y0+y1)*DEGtoRAD))
         dy = max(dy,abs(y1-y0))

!        dph  = (   x0           -     x1)*DEGtoRAD
!        dth  = sin(y0*DEGtoRAD) - sin(y1*DEGtoRAD)
!        dist = sqrt(dth**2 + dph**2)*rEarth
         dph  = (x1-x0)*cos(DEGtoRAD*(y0+y1)/2.0)
         dth  = (y1-y0)
         dist = sqrt(dth**2 + dph**2)*DEGtoRAD*rEarth
         sum = sum+dist
         maxdist = max(maxdist,dist)

       endif

   enddo
   sum = sum/sMat%n_s
   write(6,'(a,f16.3,a)') 'radius of earth at equator           :',2*pi*rEarth,'m'
   write(6,'(a,f16.3,a)') 'max distance between src & dest cells:',maxdist,'m'
   write(6,'(a,f16.3,a)') 'avg distance between src & dest cells:', sum   ,'m'
   write(6,'(a,f16.3,a)') 'max distance/radius of earth at equator:',maxdist/(2*pi*rEarth)
   write(6,'(a,f16.3,a)') 'avg distance/radius of earth at equator:',  sum  /(2*pi*rEarth)
   write(6,*) 'max dx,dy: ',dx,dy

end if

   write(6,F20) '=============================================================='

   deallocate(colsum)
   deallocate(colfsum)
   deallocate(colnsum)

END SUBROUTINE map_check

!===============================================================================

real(r8) FUNCTION map_distance(x0,y0,x1,y1)

   implicit none

   !--- arguments ---
   real(r8), intent(in) :: x0,y0,x1,y1

   !--- local ---
   real(r8) :: dph,dth,myx0,myx1

!-------------------------------------------------------------------------------
! PURPOSE:
!  o Return the distance in meters between two points on the
!    Earth with coordinates (x0,y0) and (x1,y1) in degrees.
!  o NOTE: does not take into account curvature of the grid...
!-------------------------------------------------------------------------------

   myx0 = x0
   myx1 = x1
   if (x1-x0>180.0) myx0 = x0 + 360.0
   if (x0-x1>180.0) myx1 = x1 + 360.0
   dph  = (myx1-myx0)*cos(DEGtoRAD*(y0+y1)/2.0)
   dth  = (y1-y0)
   map_distance = sqrt(dth**2 + dph**2)*DEGtoRAD*rEarth

END FUNCTION map_distance

!===============================================================================
!===============================================================================
END MODULE map_mod
!===============================================================================
