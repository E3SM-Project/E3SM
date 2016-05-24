module mkurbanparDomMod
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mkurbanparDom
!
! !DESCRIPTION:
! Make Urban Parameter data, using a dominant type approach
!
! This approach involves determining the dominant urban density class and region for each
! output grid cell, then using these indices along with lookup tables to essentially
! paint-by-number the output urban parameter fields.
!
! !REVISION HISTORY:
! Author: Bill Sacks
!
!-----------------------------------------------------------------------
! !USES:
   use shr_kind_mod, only : r8 => shr_kind_r8
   use shr_sys_mod , only : shr_sys_flush
   implicit none

   private

! !PUBLIC MEMBER FUNCTIONS:
   public :: mkurban
   public :: mkurbanpar

   ! The following could be private, but because there are associated test routines in a
   ! separate module, they need to be public
   public :: mkurban_dominant_density

! !PRIVATE DATA MEMBERS:
   ! flag to indicate nodata for index variables in output file:
   integer, parameter :: index_nodata = 0      
   character(len=*), parameter :: modname = 'mkurbanparDomMod'

   private :: index_nodata
   private :: modname
!EOP

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkurban
!
! !INTERFACE:
subroutine mkurban(ldomain, mapfname, datfname, ndiag, zero_out, urbn_o, dens_o, region_o)
!
! !DESCRIPTION:
! make percent urban, density class and region ID on the output grid
!
! Note: in contrast to mkurban in mkurbanparAvgMod, this routine does NOT handle all_urban
!
! !USES:
   use mkdomainMod , only : domain_type, domain_clean, domain_read
   use mkgridmapMod
   use mkindexmapMod, only : get_dominant_indices
   use mkurbanparCommonMod, only : mkurban_pct, mkurban_pct_diagnostics, MIN_DENS
   use mkvarctl    , only : all_urban
   use mkvarpar
   use mkncdio
!
! !ARGUMENTS:
   implicit none
   type(domain_type), intent(in) :: ldomain
   character(len=*) , intent(in) :: mapfname    ! input mapping file name
   character(len=*) , intent(in) :: datfname    ! input data file name
   integer          , intent(in) :: ndiag       ! unit number for diag out
   logical          , intent(in) :: zero_out    ! if should zero urban out
   real(r8)         , intent(out):: urbn_o(:)   ! output grid: %urban
   integer          , intent(out):: dens_o(:)   ! output grid: urban density class
   integer          , intent(out):: region_o(:) ! output grid: region ID
!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! Author: Bill Sacks
!
!
! !LOCAL VARIABLES:
!EOP
   type(domain_type) :: tdomain                ! local domain
   type(gridmap_type) :: tgridmap              ! local gridmap
   real(r8), allocatable :: urbn_by_dens_i(:,:)! input grid: percent urban in each density class, at each point
   real(r8), allocatable :: urbn_by_dens_o(:,:)! output grid: percent urban in each density class, at each point
   integer , allocatable :: region_i(:)        ! input grid: region ID
   real(r8), allocatable :: gad_i(:)           ! input grid: global area of each density class
   real(r8), allocatable :: gad_o(:)           ! output grid: global area of each density class
   real(r8), allocatable :: gar_i(:)           ! input grid: global area of each urban region ID
   real(r8), allocatable :: gar_o(:)           ! output grid: global area of each urban region ID
   real(r8) :: sum_i, sum_o                    ! sums of global areas on input & output grids
   integer  :: ni,no,ns,k                      ! indices
   integer  :: ncid,dimid,varid                ! input netCDF id's
   integer  :: dimlen                          ! netCDF dimension length
   integer  :: max_dens                        ! maximum density index
   integer  :: max_region                      ! maximum region index
   integer  :: ier                             ! error status

   character(len=32) :: subname = 'mkurban'
!-----------------------------------------------------------------------
   
   write (6,*) 'Attempting to make %urban and dominant density .....'

   ! Obtain input grid info, read local fields

   call gridmap_mapread(tgridmap, mapfname)
   call domain_read(tdomain, datfname)

   ns = tdomain%ns

   write (6,*) 'Open urban file: ', trim(datfname)
   call check_ret(nf_open(datfname, 0, ncid), subname)
   call check_ret(nf_inq_dimid (ncid, 'density_class', dimid), subname)
   call check_ret(nf_inq_dimlen (ncid, dimid, max_dens), subname)

   allocate(urbn_by_dens_i(ns, 1:max_dens), &
            urbn_by_dens_o(ldomain%ns, 1:max_dens), &
            stat=ier)
   if (ier/=0) call abort()

   call check_ret(nf_inq_varid (ncid, 'PCT_URBAN', varid), subname)
   call check_ret(nf_get_var_double (ncid, varid, urbn_by_dens_i), subname)


   ! Determine % urban by density class on the output grid
   ! Note: in some cases (e.g., zero_out=.true., or a grid cell has < MIN_DENS urban %),
   ! urbn_by_dens_o will be reset to 0 in some / all places. However, we still need the
   ! values of urbn_by_dens_o before it is zeroed, in order to compute the dominant
   ! density class in each grid cell

   do k = 1, max_dens
      ! make % urban for each density class on the output grid
      call mkurban_pct(ldomain, tdomain, tgridmap, urbn_by_dens_i(:,k), urbn_by_dens_o(:,k))
   end do


   ! Determine dominant urban density class and total % urban

   call mkurban_dominant_density(urbn_by_dens_o, index_nodata, dens_o, urbn_o)


   ! Handle special cases and too-small urban density:
   
   if (all_urban) then
      write(6,*) modname//':'//subname//' ERROR: all_urban not handled here'
      call abort()
   else if (zero_out) then
      urbn_o(:) = 0._r8
      urbn_by_dens_o(:,:) = 0._r8
   else
      do no = 1, ldomain%ns
         if (urbn_o(no) < MIN_DENS) then
            urbn_o(no) = 0._r8
            urbn_by_dens_o(no,:) = 0._r8
         end if
      end do
   end if

   ! Print diagnostics

   do k = 1, max_dens
      call mkurban_pct_diagnostics(ldomain, tdomain, tgridmap, &
                                   urbn_by_dens_i(:,k), urbn_by_dens_o(:,k), &
                                   ndiag, dens_class=k)
   end do

   write (6,*) 'Successfully made %urban and dominant density'


   write(6,*) 'Attempting to make urban region .....'

   ! Read in region field
   ! Note: we do this here, rather than with the rest of the reads above, because we
   ! expect the input urban fields to be large, so we're just reading the fields as
   ! they're needed to try to avoid unnecessary memory paging

   allocate(region_i(ns), stat=ier)
   if (ier/=0) call abort()
   call check_ret(nf_inq_varid (ncid, 'REGION_ID', varid), subname)
   call check_ret(nf_get_var_int (ncid, varid, region_i), subname)

   ! Determine max region value, and make sure it doesn't exceed bounds of the lookup tables.
   !
   ! (Note: this check assumes that region_i=1 refers to region(1), region_i=2 refers to
   ! region(2), etc. The alternative would be to use a coordinate variable associated with
   ! the region dimension of the lookup table, which could result in an arbitrary mapping
   ! between region values and the indices of the lookup table; however, this use of
   ! coordinate variables currently isn't supported by lookup_2d_netcdf [as of 2-8-12].)

   max_region = maxval(region_i)
   call check_ret(nf_inq_dimid (ncid, 'region', dimid), subname)
   call check_ret(nf_inq_dimlen (ncid, dimid, dimlen), subname)
   if (max_region > dimlen) then
      write(6,*) modname//':'//subname// &
           ' ERROR: max region value exceeds length of region dimension'
      write(6,*) 'max region value          : ', max_region
      write(6,*) 'length of region dimension: ', dimlen
      call abort()
   end if

   ! Determine dominant region for each output cell

   call get_dominant_indices(tgridmap, region_i, region_o, 1, max_region, index_nodata)

   write (6,*) 'Successfully made urban region'
   write (6,*)

   ! -----------------------------------------------------------------
   ! Error check
   ! Compare relative areas of each density class & region ID on input and output grids
   ! -----------------------------------------------------------------

   allocate(gad_i(1:max_dens),   gad_o(1:max_dens), &
            gar_i(1:max_region), gar_o(1:max_region), &
            stat=ier)
   if (ier/=0) call abort()

   gad_i(:) = 0.
   gar_i(:) = 0.
   do k = 1, max_dens
      do ni = 1,tdomain%ns
         gad_i(k) = gad_i(k) + (urbn_by_dens_i(ni,k)/100._r8)*tgridmap%area_src(ni)*tgridmap%frac_src(ni)*re**2
      end do
   end do
   do ni = 1,tdomain%ns
      k = region_i(ni)
      if (k >= 1 .and. k <= max_region) then
         gar_i(k) = gar_i(k) + tgridmap%area_src(ni)*tgridmap%frac_src(ni)*re**2
      end if
   end do

   gad_o(:) = 0.
   gar_o(:) = 0.
   do no = 1,ldomain%ns
      k = dens_o(no)
      if (k >= 1 .and. k <= max_dens) then
         gad_o(k) = gad_o(k) + tgridmap%area_dst(no)*tgridmap%frac_dst(no)*re**2
      end if
      k = region_o(no)
      if (k >= 1 .and. k <= max_region) then
         gar_o(k) = gar_o(k) + tgridmap%area_dst(no)*tgridmap%frac_dst(no)*re**2
      end if      
   end do
  
   ! relative area comparison
   ! Note: we compare relative areas rather than absolute areas because we expect big
   ! differences in the absolute areas between input & output.
   ! Also note that the relative areas are as a proportion of the area with valid urban data
   
   write (ndiag,*)
   write (ndiag,'(1x,70a1)') ('=',k=1,70)
   write (ndiag,*) 'Urban Density Class Output'
   write (ndiag,'(1x,70a1)') ('=',k=1,70)

   write (ndiag,*)
   write (ndiag,'(1x,70a1)') ('.',k=1,70)
   write (ndiag,1001)
1001 format (1x,'density class input grid area output grid area',/ &
             1x,'                      percent          percent')
   write (ndiag,'(1x,70a1)') ('.',k=1,70)
   write (ndiag,*)

   sum_i = sum(gad_i)
   sum_o = sum(gad_o)
   do k = 1, max_dens
      write (ndiag,1002) k,gad_i(k)/sum_i*100,gad_o(k)/sum_o*100
1002  format (1x,i13,f15.3,'%',f16.3,'%')
   end do

   write (ndiag,*)
   write (ndiag,'(1x,70a1)') ('=',k=1,70)
   write (ndiag,*) 'Urban Region ID Output'
   write (ndiag,'(1x,70a1)') ('=',k=1,70)

   write (ndiag,*)
   write (ndiag,'(1x,70a1)') ('.',k=1,70)
   write (ndiag,1003)
1003 format (1x,'region ID input grid area output grid area',/ &
             1x,'                  percent          percent')
   write (ndiag,'(1x,70a1)') ('.',k=1,70)
   write (ndiag,*)

   sum_i = sum(gar_i)
   sum_o = sum(gar_o)
   do k = 1, max_region
      write (ndiag,1004) k,gar_i(k)/sum_i*100,gar_o(k)/sum_o*100
1004  format (1x,i9,f15.3,'%',f16.3,'%')
   end do


   ! Deallocate dynamic memory & other clean up

   call check_ret(nf_close(ncid), subname)
   call domain_clean(tdomain)
   call gridmap_clean(tgridmap)
   deallocate (urbn_by_dens_i, urbn_by_dens_o, region_i)
   deallocate(gad_i, gad_o, gar_i, gar_o)
  
end subroutine mkurban
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkurbanpar
!
! !INTERFACE:
subroutine mkurbanpar(datfname, ncido, dens_o, region_o)
!
! !DESCRIPTION:
! Make Urban Parameter data
!
! !USES:
   use mkdomainMod  , only : domain_type, domain_clean, domain_read
   use mkindexmapMod, only : dim_slice_type, lookup_2d_netcdf
   use mkvarpar
   use mkncdio
!
! !ARGUMENTS:
   implicit none
   character(len=*)  , intent(in) :: datfname    ! input data file name
   integer           , intent(in) :: ncido       ! output netcdf file id
   integer           , intent(in) :: dens_o(:)   ! output grid: urban density class
   integer           , intent(in) :: region_o(:) ! output grid: region ID
!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! Author: Bill Sacks
!
!
! !LOCAL VARIABLES:
!EOP
   type param
      character(len=nf_max_name) :: name         ! name in input & output files
      logical                    :: invalid_okay ! are NA values allowed in the input table?
   end type param

   real(r8), allocatable :: data_scalar_o(:)   ! output array for parameters with no extra dimensions
   real(r8), allocatable :: data_rad_o(:,:,:)  ! output array for parameters dimensioned by numrad & numsolar
   real(r8), allocatable :: data_levurb_o(:,:) ! output array for parameters dimensioned by nlevurb
   integer  :: nlevurb_i                       ! input  grid: number of urban vertical levels
   integer  :: numsolar_i                      ! input  grid: number of solar type (DIR/DIF)
   integer  :: numrad_i                        ! input  grid: number of solar bands (VIS/NIR)
   integer  :: m,n,no,ns_o,p                   ! indices
   integer  :: ncidi,dimid,varid               ! netCDF id's
   integer  :: ier                             ! error status
   
   ! information on extra dimensions for lookup tables greater than 2-d:
   type(dim_slice_type), allocatable :: extra_dims(:)  

   ! To add a new urban parameter, simply add an element to one of the below lists
   ! (params_scalar, params_rad or params_levurb)
   
   ! Urban parameters with no extra dimensions
   type(param), parameter :: params_scalar(14) = &
        (/ param('CANYON_HWR', .false.), &
           param('EM_IMPROAD', .false.), &
           param('EM_PERROAD', .false.), &
           param('EM_ROOF', .false.), &
           param('EM_WALL', .false.), &
           param('HT_ROOF', .false.), &
           param('THICK_ROOF', .false.), &
           param('THICK_WALL', .false.), &
           param('T_BUILDING_MAX', .false.), &
           param('T_BUILDING_MIN', .false.), &
           param('WIND_HGT_CANYON', .false.), &
           param('WTLUNIT_ROOF', .false.), &
           param('WTROAD_PERV', .false.), &

           ! Note that NLEV_IMPROAD is written as an integer, meaning that type conversion occurs
           ! by truncation. Thus we expect the values in the NLEV_IMPROAD lookup table to be exact;
           ! e.g., if a value were 1.99999 rather than 2.0000, it would be written as 1 instead of 2
           param('NLEV_IMPROAD', .false.) /)

   ! Urban parameters dimensioned by numrad & numsolar
   type(param), parameter :: params_rad(4) = &
        (/ param('ALB_IMPROAD', .false.), &
           param('ALB_PERROAD', .false.), &
           param('ALB_ROOF', .false.), &
           param('ALB_WALL', .false.) /)
   
   ! Urban parameters dimensioned by nlevurb
   type(param), parameter :: params_levurb(6) = &
        (/ param('TK_ROOF', .false.), &
           param('TK_WALL', .false.), &
           param('CV_ROOF', .false.), &
           param('CV_WALL', .false.), &

           ! Impervious road thermal conductivity and heat capacity have varying levels of
           ! data. Thus, we expect to find some missing values in the lookup table -- we
           ! do not want to treat that as an error -- thus, we set invalid_okay=.true.
           param('CV_IMPROAD', .true.), &
           param('TK_IMPROAD', .true.) /)


   character(len= 32) :: subname = 'mkurbanpar'
!-----------------------------------------------------------------------

   write (6,*) 'Attempting to make Urban Parameters .....'
   call shr_sys_flush(6)

   ! Determine & error-check array sizes
   ns_o = size(dens_o)
   if (ns_o /= size(region_o)) then
      write(6,*) modname//':'//subname//' ERROR: array size mismatch'
      write(6,*) 'size(dens_o)   = ', size(dens_o)
      write(6,*) 'size(region_o) = ', size(region_o)
      call abort()
   end if

   ! Read dimensions from input file

   write (6,*) 'Open urban parameter file: ', trim(datfname)
   call check_ret(nf_open(datfname, 0, ncidi), subname)
   call check_ret(nf_inq_dimid(ncidi, 'nlevurb', dimid), subname)
   call check_ret(nf_inq_dimlen(ncidi, dimid, nlevurb_i), subname)
   call check_ret(nf_inq_dimid(ncidi, 'numsolar', dimid), subname)
   call check_ret(nf_inq_dimlen(ncidi, dimid, numsolar_i), subname)
   call check_ret(nf_inq_dimid(ncidi, 'numrad', dimid), subname)
   call check_ret(nf_inq_dimlen(ncidi, dimid, numrad_i), subname)

   if (nlevurb_i /= nlevurb) then
      write(6,*)'MKURBANPAR: parameter nlevurb= ',nlevurb, &
           'does not equal input dataset nlevurb= ',nlevurb_i
      stop
   endif
   if (numsolar_i /= numsolar) then
      write(6,*)'MKURBANPAR: parameter numsolar= ',numsolar, &
           'does not equal input dataset numsolar= ',numsolar_i
      stop
   endif
   if (numrad_i /= numrad) then
      write(6,*)'MKURBANPAR: parameter numrad= ',numrad, &
           'does not equal input dataset numrad= ',numrad_i
      stop
   endif

   ! Handle urban parameters with no extra dimensions

   allocate(data_scalar_o(ns_o), stat=ier)
   if (ier /= 0) then
      write(6,*)'mkurbanpar allocation error'; call abort()
   end if
   
   do p = 1, size(params_scalar)
      call lookup_and_check_err(params_scalar(p)%name, data_scalar_o, 0, &
                                invalid_okay=params_scalar(p)%invalid_okay)

      call check_ret(nf_inq_varid(ncido, params_scalar(p)%name, varid), subname)
      ! In the following, note that type conversion occurs if we're writing to a variable of type
      ! other than double; e.g., for an integer, conversion occurs by truncation!
      call check_ret(nf_put_var_double(ncido, varid, data_scalar_o), subname)
   end do

   deallocate(data_scalar_o)
      
   ! Handle urban parameters dimensioned by numrad & numsolar

   allocate(data_rad_o(ns_o, numrad, numsolar), stat=ier)
   if (ier /= 0) then
      write(6,*)'mkurbanpar allocation error'; call abort()
   end if

   allocate(extra_dims(2))
   extra_dims(1)%name = 'numrad'
   extra_dims(2)%name = 'numsolar'

   do p = 1, size(params_rad)
      do m = 1,numsolar_i
         extra_dims(2)%val = m
         do n = 1,numrad_i
            extra_dims(1)%val = n

            call lookup_and_check_err(params_rad(p)%name, data_rad_o(:,n,m), 2, extra_dims, &
                                      invalid_okay=params_rad(p)%invalid_okay)
         end do
      end do

      call check_ret(nf_inq_varid(ncido, params_rad(p)%name, varid), subname)
      ! In the following, note that type conversion occurs if we're writing to a variable of type
      ! other than double; e.g., for an integer, conversion occurs by truncation!
      call check_ret(nf_put_var_double(ncido, varid, data_rad_o), subname)
   end do

   deallocate(data_rad_o)
   deallocate(extra_dims)

   ! Handle urban parameters dimensioned by nlevurb

   allocate(data_levurb_o(ns_o, nlevurb), stat=ier)
   if (ier /= 0) then
      write(6,*)'mkurbanpar allocation error'; call abort()
   end if   

   allocate(extra_dims(1))
   extra_dims(1)%name = 'nlevurb'

   do p = 1, size(params_levurb)
      do n = 1,nlevurb
         extra_dims(1)%val = n

         call lookup_and_check_err(params_levurb(p)%name, data_levurb_o(:,n), 1, extra_dims, &
                                   invalid_okay=params_levurb(p)%invalid_okay)
      end do

      call check_ret(nf_inq_varid(ncido, params_levurb(p)%name, varid), subname)
      ! In the following, note that type conversion occurs if we're writing to a variable of type
      ! other than double; e.g., for an integer, conversion occurs by truncation!
      call check_ret(nf_put_var_double(ncido, varid, data_levurb_o), subname)      
   end do

   deallocate(data_levurb_o)
   deallocate(extra_dims)
   

   call check_ret(nf_close(ncidi), subname)
   call check_ret(nf_sync(ncido), subname)

   write (6,*) 'Successfully made Urban Parameters'
   write (6,*)
   call shr_sys_flush(6)

contains
!------------------------------------------------------------------------------
   subroutine lookup_and_check_err(varname, data, n_extra_dims, extra_dims, invalid_okay)
   ! Wrapper to lookup_2d_netcdf: calls that routine with the appropriate arguments,
   ! then checks the error code, aborting if there were errors
   !
   ! Note: inherits a number of variables from the parent routine

      use mkindexmapMod, only : lookup_2d_netcdf

      implicit none
      character(len=*), intent(in) :: varname      ! name of lookup table
      real(r8)        , intent(out):: data(:)      ! output from lookup_2d_netcdf
      integer         , intent(in) :: n_extra_dims ! number of extra dimensions in the lookup table

      ! slice to use if lookup table variable has more than 2 dimensions:
      type(dim_slice_type), intent(in), optional :: extra_dims(:)

      ! if present and true, then we won't abort due to finding _FillValue in the lookup
      ! table -- instead, those places will just have fill_val in the output data
      logical, intent(in), optional :: invalid_okay

      ! Local variables:

      integer :: ierr  ! error return code

      ! value to put where we have no data in output variables
      real(r8), parameter :: fill_val = 0._r8


      call lookup_2d_netcdf(ncidi, varname, .true., &
                            'density_class', 'region', n_extra_dims, &
                            dens_o, region_o, fill_val, data, ierr, &
                            extra_dims=extra_dims, nodata=index_nodata, &
                            invalid_okay=invalid_okay)

      if (ierr /= 0) then
         write(6,*) modname//':'//subname//' ERROR in lookup_2d_netcdf for ', &
              trim(varname), ':', ierr
         call abort()
      end if
   end subroutine lookup_and_check_err

end subroutine mkurbanpar
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkurban_dominant_density
!
! !INTERFACE:
subroutine mkurban_dominant_density(urbn_by_dens_o, nodata, dens_o, urbn_o)
!
! !DESCRIPTION:
! Creates urban density class and total % urban on the output grid
!
! !USES:
   use mkindexmapMod, only : which_max
!
! !ARGUMENTS:
   implicit none
   real(r8), intent(in) :: urbn_by_dens_o(:,:) ! % urban in each density class, output grid
                                               ! (dimensions: (gridmap%nb, ndens))
   integer , intent(in) :: nodata              ! flag to indicate nodata in output arrays
   integer , intent(out):: dens_o(:)           ! urban density class, output grid
   real(r8), intent(out):: urbn_o(:)           ! total % urban, output grid
! 
! !REVISION HISTORY:
! Author: Bill Sacks
!
!
! !LOCAL VARIABLES:
!EOP
   integer  :: n
   integer  :: no
   real(r8) :: maxval
   integer  :: maxindex
   
   character(len=*), parameter :: subName = "mkurban_dominant_density"
!------------------------------------------------------------------------------

   n = size(urbn_by_dens_o, 1)

   ! Error-check inputs

   if (size(dens_o) /= n .or. size(urbn_o) /= n) then
      write(6,*) subName//' ERROR: incorrect array sizes'
      write(6,*) 'n            = ', n
      write(6,*) 'size(dens_o) = ', size(dens_o)
      write(6,*) 'size(urbn_o) = ', size(urbn_o)
      call abort()
   end if


   do no = 1, size(dens_o)
      urbn_o(no) = sum(urbn_by_dens_o(no,:))

      ! Determine dominant density class for each output cell
      ! Note: if all urban density classes have 0 area, then the output value will be nodata
      call which_max(urbn_by_dens_o(no,:), maxval, maxindex)
      if (maxval > 0.) then
         dens_o(no) = maxindex
      else
         dens_o(no) = nodata
      end if
   end do

end subroutine mkurban_dominant_density
!------------------------------------------------------------------------------

end module mkurbanparDomMod
