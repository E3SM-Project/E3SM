module mkurbanparMod
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mkurbanpar
!
! !DESCRIPTION:
! Make Urban Parameter data
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
   public :: mkurbanInit
   public :: mkurban
   public :: mkurbanpar

   ! The following could be private, but because there are associated test routines in a
   ! separate module, it needs to be public
   public :: normalize_urbn_by_tot

! !PUBLIC DATA MEMBERS:
   integer :: numurbl           ! number of urban classes

   public :: numurbl

! !PRIVATE DATA MEMBERS:
   ! flag to indicate nodata for index variables in output file:
   integer, parameter :: index_nodata = 0      
   character(len=*), parameter :: modname = 'mkurbanparMod'

   private :: index_nodata
   private :: modname

!EOP

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkurbanInit
!
! !INTERFACE:
subroutine mkurbanInit(datfname)
!
! !DESCRIPTION:
! Initialize variables needed for urban
!
! !USES:
   use mkncdio
!
! !ARGUMENTS:
   implicit none
   character(len=*), intent(in) :: datfname  ! input data file name (same as file used in mkurban)
!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! Author: Bill Sacks
!
!
! !LOCAL VARIABLES:
   integer  :: ncid,dimid                       ! input netCDF id's

   character(len=*), parameter :: subname = 'mkurbanInit'
!EOP
!-----------------------------------------------------------------------

   ! Set numurbl
   call check_ret(nf_open(datfname, 0, ncid), subname)
   call check_ret(nf_inq_dimid (ncid, 'density_class', dimid), subname)
   call check_ret(nf_inq_dimlen (ncid, dimid, numurbl), subname)
   call check_ret(nf_close(ncid), subname)

end subroutine mkurbanInit
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkurban
!
! !INTERFACE:
subroutine mkurban(ldomain, mapfname, datfname, ndiag, zero_out, &
                   urbn_o, urbn_classes_o, region_o)
!
! !DESCRIPTION:
! make total percent urban, breakdown into urban classes, and region ID on the output grid
!
! urbn_classes_o(n, i) gives the percent of the urban area in grid cell n that is in class #i.
! This is normalized so that sum(urbn_classes_o(n,:)) = 100 for all n, even for grid
! cells where urbn_o(n) = 0 (in the case where urbn_o(n) = 0, we come up with an
! arbitrary assignment of urban into the different classes).
!
! See comments under the normalize_urbn_by_tot subroutine for how urbn_classes_o is
! determined when the total % urban is 0, according to the input data. Note that this
! also applies when all_urban=.true., for points that have 0 urban according to the input
! data.
!
! !USES:
   use mkdomainMod , only : domain_type, domain_clean, domain_read
   use mkgridmapMod
   use mkindexmapMod, only : get_dominant_indices
   use mkurbanparCommonMod, only : mkurban_pct, mkurban_pct_diagnostics, MIN_DENS
   use mkutilsMod  , only : normalize_classes_by_gcell
   use mkvarctl    , only : all_urban
   use mkvarpar
   use mkncdio
!
! !ARGUMENTS:
   implicit none
   type(domain_type), intent(in) :: ldomain
   character(len=*) , intent(in) :: mapfname            ! input mapping file name
   character(len=*) , intent(in) :: datfname            ! input data file name
   integer          , intent(in) :: ndiag               ! unit number for diag out
   logical          , intent(in) :: zero_out            ! if should zero urban out
   real(r8)         , intent(out):: urbn_o(:)           ! output grid: total % urban
   real(r8)         , intent(out):: urbn_classes_o(:,:) ! output grid: breakdown of total urban into each class
                                                        ! (dimensions: (ldomain%ns, numurbl))
   integer          , intent(out):: region_o(:)         ! output grid: region ID
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
   type(domain_type) :: tdomain                       ! local domain
   type(gridmap_type) :: tgridmap                     ! local gridmap
   real(r8), allocatable :: urbn_classes_gcell_i(:,:) ! input grid: percent urban in each density class
                                                      ! (% of total grid cell area)
   real(r8), allocatable :: urbn_classes_gcell_o(:,:) ! output grid: percent urban in each density class
                                                      ! (% of total grid cell area)
   integer , allocatable :: region_i(:)               ! input grid: region ID
   real(r8), allocatable :: gar_i(:)                  ! input grid: global area of each urban region ID
   real(r8), allocatable :: gar_o(:)                  ! output grid: global area of each urban region ID
   integer  :: ni,no,ns,k                             ! indices
   integer  :: ncid,dimid,varid                       ! input netCDF id's
   integer  :: dimlen                                 ! netCDF dimension length
   integer  :: max_region                             ! maximum region index
   integer  :: ier                                    ! error status

   character(len=*), parameter :: subname = 'mkurban'
!-----------------------------------------------------------------------
   
   write (6,*) 'Attempting to make %urban .....'

   ! Obtain input grid info, read local fields

   call gridmap_mapread(tgridmap, mapfname)
   call domain_read(tdomain, datfname)

   ns = tdomain%ns

   allocate(urbn_classes_gcell_i(ns, numurbl), &
            urbn_classes_gcell_o(ldomain%ns, numurbl), &
            stat=ier)
   if (ier/=0) call abort()

   write (6,*) 'Open urban file: ', trim(datfname)
   call check_ret(nf_open(datfname, 0, ncid), subname)
   call check_ret(nf_inq_varid (ncid, 'PCT_URBAN', varid), subname)
   call check_ret(nf_get_var_double (ncid, varid, urbn_classes_gcell_i), subname)

   ! Determine % urban by density class on the output grid
   do k = 1, numurbl
      call mkurban_pct(ldomain, tdomain, tgridmap, urbn_classes_gcell_i(:,k), urbn_classes_gcell_o(:,k))
   end do

   ! Determine total % urban
   do no = 1, ldomain%ns
      urbn_o(no) = sum(urbn_classes_gcell_o(no,:))
   end do

   call normalize_urbn_by_tot(urbn_classes_gcell_o, urbn_o, urbn_classes_o)

   ! Handle special cases

   ! Note that, for all these adjustments of total urban %, we do not change anything
   ! about the breakdown into the different urban classes. In particular: when urbn_o is
   ! set to 0 for a point, the breakdown into the different urban classes is maintained
   ! as it was before.
   if (all_urban) then
      urbn_o(:) = 100._r8
   else if (zero_out) then
      urbn_o(:) = 0._r8
   else
      ! Set points to 0% if they fall below a given threshold
      do no = 1, ldomain%ns
         if (urbn_o(no) < MIN_DENS) then
            urbn_o(no) = 0._r8
         end if
      end do
   end if

   ! Print diagnostics
   ! First, recompute urbn_classes_gcell_o, based on any changes we have made to urbn_o
   ! while handling special cases
   call normalize_classes_by_gcell(urbn_classes_o, urbn_o, urbn_classes_gcell_o)
   do k = 1, numurbl
      call mkurban_pct_diagnostics(ldomain, tdomain, tgridmap, &
           urbn_classes_gcell_i(:,k), urbn_classes_gcell_o(:,k), &
           ndiag, dens_class=k)
   end do

   write (6,*) 'Successfully made %urban'


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
   ! Compare areas of each region ID on input and output grids
   ! -----------------------------------------------------------------

   allocate(gar_i(max_region), gar_o(max_region), stat=ier)
   if (ier/=0) call abort()

   gar_i(:) = 0.
   do ni = 1,tdomain%ns
      k = region_i(ni)
      if (k >= 1 .and. k <= max_region) then
         gar_i(k) = gar_i(k) + tgridmap%area_src(ni)*tgridmap%frac_src(ni)*re**2
      end if
   end do

   gar_o(:) = 0.
   do no = 1,ldomain%ns
      k = region_o(no)
      if (k >= 1 .and. k <= max_region) then
         gar_o(k) = gar_o(k) + tgridmap%area_dst(no)*tgridmap%frac_dst(no)*re**2
      end if      
   end do
  
   write (ndiag,*)
   write (ndiag,'(1x,70a1)') ('=',k=1,70)
   write (ndiag,*) 'Urban Region ID Output'
   write (ndiag,'(1x,70a1)') ('=',k=1,70)

   write (ndiag,*)
   write (ndiag,'(1x,70a1)') ('.',k=1,70)
   write (ndiag,1003)
1003 format (1x,'region ID  input grid area  output grid area',/ &
             1x,'               10**6 km**2       10**6 km**2')
   write (ndiag,'(1x,70a1)') ('.',k=1,70)
   write (ndiag,*)

   do k = 1, max_region
      write (ndiag,1004) k,gar_i(k)*1.e-06,gar_o(k)*1.e-06
1004  format (1x,i9,f17.3,f18.3)
   end do


   ! Deallocate dynamic memory & other clean up

   call check_ret(nf_close(ncid), subname)
   call domain_clean(tdomain)
   call gridmap_clean(tgridmap)
   deallocate (urbn_classes_gcell_i, urbn_classes_gcell_o, region_i)
   deallocate(gar_i, gar_o)
  
end subroutine mkurban
!-----------------------------------------------------------------------

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: normalize_urbn_by_tot
!
! !INTERFACE:
subroutine normalize_urbn_by_tot(classes_pct_gcell, sums, classes_pct_tot)
!
! !DESCRIPTION:
! Normalizes urban class areas to produce % cover of each class, as % of total urban area
!
! Specifically: Given (1) an array specifying the % cover of each urban class, as a % of
! the total grid cell area ('classes_pct_gcell'), and (2) a vector giving the total urban
! area in each grid cell, expressed as % of the grid cell area: Returns an array
! ('classes_pct_tot') of the same dimensionality as classes_pct_gcell, where the values
! now give % cover of each class as a % of the total urban area.
!
! Assumes that sums(n) = sum(classes_pct_gcell(n,:))
!
! When sums(n) = 0, the creation of classes_pct_tot(n,:) is ambiguous. Here we use the
! rule that all area is assigned to the medium-density class, defined by parameter MD.
!
! The returned array satisfies sum(classes_pct_tot(n,:))==100 for all n (within rounding error)
!
! !USES:
!
! !ARGUMENTS:
   implicit none
   real(r8), intent(in) :: classes_pct_gcell(:,:) ! % cover of classes as % of grid cell
   real(r8), intent(in) :: sums(:)                ! totals, as % of grid cell
   real(r8), intent(out):: classes_pct_tot(:,:)   ! % cover of classes as % of total
!
! !REVISION HISTORY:
! Author: Bill Sacks
!
!
! !LOCAL VARIABLES:
!EOP
   integer  :: n         ! index
   integer  :: n_max     ! number of points 
   integer  :: nclasses  ! number of classes
   real(r8) :: suma      ! sum for error check
   
   ! index of medium-density class, which is where we assign urban areas when the total
   ! urban area is 0
   integer, parameter :: MD = 3

   ! relative error tolerance for error check
   real(r8), parameter :: relerr = 1.e-10_r8
   
   character(len=*), parameter :: subname = 'normalize_urbn_by_tot'
!-----------------------------------------------------------------------

   ! Error-check inputs

   n_max = size(sums)
   if (size(classes_pct_tot, 1)   /= n_max .or. &
       size(classes_pct_gcell, 1) /= n_max) then
      write(6,*) subname//' ERROR: array size mismatch'
      write(6,*) 'size(sums)                 = ', n_max
      write(6,*) 'size(classes_pct_tot, 1)   = ', size(classes_pct_tot, 1)
      write(6,*) 'size(classes_pct_gcell, 1) = ', size(classes_pct_gcell, 1)
      call abort()
   end if

   if (size(classes_pct_tot, 2) /= size(classes_pct_gcell, 2)) then
      write(6,*) subname//' ERROR: array size mismatch'
      write(6,*) 'size(classes_pct_tot, 2)   = ', size(classes_pct_tot, 2)
      write(6,*) 'size(classes_pct_gcell, 2) = ', size(classes_pct_gcell, 2)
      call abort()
   end if

   nclasses = size(classes_pct_gcell, 2)
   if (MD > nclasses) then
      write(6,*) subname//' ERROR: MD exceeds nclasses'
      write(6,*) 'MD       = ', MD
      write(6,*) 'nclasses = ', nclasses
      call abort()
   end if

   ! Do the work

   do n = 1, n_max
      if (sums(n) > 0._r8) then
         classes_pct_tot(n,:) = classes_pct_gcell(n,:)/sums(n) * 100._r8
      else
         ! Creation of classes_pct_tot is ambiguous. Apply the rule that all area is
         ! assigned to the medium-density class.
         classes_pct_tot(n,:)  =   0._r8
         classes_pct_tot(n,MD) = 100._r8
      end if
   end do

   ! Error-check output: Make sure sum(classes_pct_tot(n,:)) = 100 for all n
   
   do n = 1, n_max
      suma = sum(classes_pct_tot(n,:))
      if (abs(suma/100._r8 - 1._r8) > relerr) then
         write(6,*) subname//' ERROR: sum does not equal 100 at point ', n
         write(6,*) 'suma = ', suma
         call abort()
      end if
   end do

end subroutine normalize_urbn_by_tot
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkurbanpar
!
! !INTERFACE:
subroutine mkurbanpar(datfname, ncido, region_o, urbn_classes_gcell_o)
!
! !DESCRIPTION:
! Make Urban Parameter data
!
! Note that, in a grid cell with region_o==r, parameter values are filled from region r
! for ALL density classes. Thus, the parameter variables have a numurbl dimension along
! with their other dimensions.
!
! Note that we will have a 'nodata' value (given by the fill_val value associated with
! each parameter) wherever (1) we have a nodata value for region_o, or (2) the parameter
! has nodata for the given region/density combination in the input lookup table.
!
! !USES:
   use mkdomainMod  , only : domain_type, domain_clean, domain_read
   use mkindexmapMod, only : dim_slice_type, lookup_2d_netcdf
   use mkvarpar
   use mkncdio
!
! !ARGUMENTS:
   implicit none
   character(len=*)  , intent(in) :: datfname                  ! input data file name
   integer           , intent(in) :: ncido                     ! output netcdf file id
   integer           , intent(in) :: region_o(:)               ! output grid: region ID (length: ns_o)
   real(r8)          , intent(in) :: urbn_classes_gcell_o(:,:) ! output grid: percent urban in each density class
                                                               ! (% of total grid cell area) (dimensions: ns_o, numurbl)

! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! Author: Bill Sacks
!
!
! !LOCAL VARIABLES:
!EOP
   ! Type to store information about each urban parameter
   type param
      character(len=32) :: name          ! name in input & output files
      real(r8)          :: fill_val      ! value to put where we have no data in output
      logical           :: check_invalid ! should we check whether there are any invalid data in the output?
   end type param

   real(r8), allocatable :: data_scalar_o(:,:)   ! output array for parameters with no extra dimensions
   real(r8), allocatable :: data_rad_o(:,:,:,:)  ! output array for parameters dimensioned by numrad & numsolar
   real(r8), allocatable :: data_levurb_o(:,:,:) ! output array for parameters dimensioned by nlevurb
   integer , allocatable :: unity_dens_o(:,:)    ! artificial density indices
   integer  :: nlevurb_i                         ! input  grid: number of urban vertical levels
   integer  :: numsolar_i                        ! input  grid: number of solar type (DIR/DIF)
   integer  :: numrad_i                          ! input  grid: number of solar bands (VIS/NIR)
   integer  :: m,n,no,ns_o,p,k                   ! indices
   integer  :: ncidi,dimid,varid                 ! netCDF id's
   integer  :: ier                               ! error status
   character(len=nf_max_name) :: varname         ! variable name
   
   ! information on extra dimensions for lookup tables greater than 2-d:
   type(dim_slice_type), allocatable :: extra_dims(:)  
   
   ! suffix for variables dimensioned by numsolar, for each value of numsolar:
   character(len=8), parameter :: solar_suffix(numsolar) = (/'_DIR', '_DIF'/)

   ! value to put where we have no data in output variables, for real-valued parameters
   real(r8), parameter :: fill_val_real = 0._r8

   ! To add a new urban parameter, simply add an element to one of the below lists
   ! (params_scalar, params_rad or params_levurb)
   
   ! Urban parameters with no extra dimensions
   type(param), parameter :: params_scalar(14) = &
        (/ param('CANYON_HWR', fill_val_real, .true.), &
           param('EM_IMPROAD', fill_val_real, .true.), &
           param('EM_PERROAD', fill_val_real, .true.), &
           param('EM_ROOF', fill_val_real, .true.), &
           param('EM_WALL', fill_val_real, .true.), &
           param('HT_ROOF', fill_val_real, .true.), &
           param('THICK_ROOF', fill_val_real, .true.), &
           param('THICK_WALL', fill_val_real, .true.), &
           param('T_BUILDING_MAX', fill_val_real, .true.), &
           param('T_BUILDING_MIN', fill_val_real, .true.), &
           param('WIND_HGT_CANYON', fill_val_real, .true.), &
           param('WTLUNIT_ROOF', fill_val_real, .true.), &
           param('WTROAD_PERV', fill_val_real, .true.), &

           ! Note that NLEV_IMPROAD is written as an integer, meaning that type conversion occurs
           ! by truncation. Thus we expect the values in the NLEV_IMPROAD lookup table to be exact;
           ! e.g., if a value were 1.99999 rather than 2.0000, it would be written as 1 instead of 2
           ! Also note: we use fill_val=-1 rather than 0, because 0 appears in the lookup table
           param('NLEV_IMPROAD', -1, .true.) /)

   ! Urban parameters dimensioned by numrad & numsolar
   type(param), parameter :: params_rad(4) = &
        (/ param('ALB_IMPROAD', fill_val_real, .true.), &
           param('ALB_PERROAD', fill_val_real, .true.), &
           param('ALB_ROOF', fill_val_real, .true.), &
           param('ALB_WALL', fill_val_real, .true.) /)
   
   ! Urban parameters dimensioned by nlevurb
   type(param), parameter :: params_levurb(6) = &
        (/ param('TK_ROOF', fill_val_real, .true.), &
           param('TK_WALL', fill_val_real, .true.), &
           param('CV_ROOF', fill_val_real, .true.), &
           param('CV_WALL', fill_val_real, .true.), &

           ! Impervious road thermal conductivity and heat capacity have varying levels of
           ! data. Thus, we expect to find some missing values in the lookup table -- we
           ! do not want to treat that as an error -- thus, we set check_invalid=.false.
           param('CV_IMPROAD', fill_val_real, .false.), &
           param('TK_IMPROAD', fill_val_real, .false.) /)


   character(len=*), parameter :: subname = 'mkurbanpar'
!-----------------------------------------------------------------------

   write (6,*) 'Attempting to make Urban Parameters .....'
   call shr_sys_flush(6)

   ! Determine & error-check array sizes
   ns_o = size(region_o)
   if (size(urbn_classes_gcell_o, 1) /= ns_o) then
      write(6,*) modname//':'//subname//' ERROR: array size mismatch'
      write(6,*) 'size(region_o) = ', size(region_o)
      write(6,*) 'size(urbn_classes_gcell_o, 1) = ', size(urbn_classes_gcell_o, 1)
      call abort()
   end if
   if (size(urbn_classes_gcell_o, 2) /= numurbl) then
      write(6,*) modname//':'//subname//' ERROR: array size mismatch'
      write(6,*) 'size(urbn_classes_gcell_o, 2) = ', size(urbn_classes_gcell_o, 2)
      write(6,*) 'numurbl = ', numurbl
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

   ! Create an array that will hold the density indices
   ! In a given grid cell, we output parameter values for all density classes, for the
   ! region of that grid cell. In order to do this while still using the lookup_2d
   ! routine, we create a dummy unity_dens_o array that contains the density values
   ! passed to the lookup routine.

   allocate(unity_dens_o(ns_o, numurbl))
   do k = 1, numurbl
      unity_dens_o(:,k) = k
   end do

   ! Handle urban parameters with no extra dimensions

   allocate(data_scalar_o(ns_o, numurbl), stat=ier)
   if (ier /= 0) then
      write(6,*)'mkurbanpar allocation error'; call abort()
   end if
   
   do p = 1, size(params_scalar)
      call lookup_and_check_err(params_scalar(p)%name, params_scalar(p)%fill_val, &
                                params_scalar(p)%check_invalid, data_scalar_o, 0)

      call check_ret(nf_inq_varid(ncido, params_scalar(p)%name, varid), subname)
      ! In the following, note that type conversion occurs if we're writing to a variable of type
      ! other than double; e.g., for an integer, conversion occurs by truncation!
      call check_ret(nf_put_var_double(ncido, varid, data_scalar_o), subname)
   end do

   deallocate(data_scalar_o)
      
   ! Handle urban parameters dimensioned by numrad & numsolar

   allocate(data_rad_o(ns_o, numurbl, numrad, numsolar), stat=ier)
   if (ier /= 0) then
      write(6,*)'mkurbanpar allocation error'; call abort()
   end if

   allocate(extra_dims(2))
   extra_dims(1)%name = 'numrad'
   extra_dims(2)%name = 'numsolar'

   do p = 1, size(params_rad)
      do m = 1,numsolar
         extra_dims(2)%val = m
         do n = 1,numrad
            extra_dims(1)%val = n

            call lookup_and_check_err(params_rad(p)%name, params_rad(p)%fill_val, &
                                      params_rad(p)%check_invalid, data_rad_o(:,:,n,m), &
                                      2, extra_dims)
         end do
      end do

      ! Special handling of numsolar: rather than outputting variables with a numsolar
      ! dimension, we output separate variables for each value of numsolar
      do m = 1,numsolar
         if (len_trim(params_rad(p)%name) + len_trim(solar_suffix(m)) > len(varname)) then
            write(6,*) 'variable name exceeds length of varname'
            write(6,*) trim(params_rad(p)%name)//trim(solar_suffix(m))
            call abort()
         end if
         varname = trim(params_rad(p)%name)//trim(solar_suffix(m))
         call check_ret(nf_inq_varid(ncido, varname, varid), subname)
         ! In the following, note that type conversion occurs if we're writing to a variable of type
         ! other than double; e.g., for an integer, conversion occurs by truncation!
         call check_ret(nf_put_var_double(ncido, varid, data_rad_o(:,:,:,m)), subname)
      end do
   end do

   deallocate(data_rad_o)
   deallocate(extra_dims)

   ! Handle urban parameters dimensioned by nlevurb

   allocate(data_levurb_o(ns_o, numurbl, nlevurb), stat=ier)
   if (ier /= 0) then
      write(6,*)'mkurbanpar allocation error'; call abort()
   end if   

   allocate(extra_dims(1))
   extra_dims(1)%name = 'nlevurb'

   do p = 1, size(params_levurb)
      do n = 1,nlevurb
         extra_dims(1)%val = n

         call lookup_and_check_err(params_levurb(p)%name, params_levurb(p)%fill_val, &
                                   params_levurb(p)%check_invalid, data_levurb_o(:,:,n), &
                                   1, extra_dims)
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

   deallocate(unity_dens_o)

contains
!------------------------------------------------------------------------------
   subroutine lookup_and_check_err(varname, fill_val, check_invalid, data, n_extra_dims, extra_dims)
   ! Wrapper to lookup_2d_netcdf: Loops over each density class, calling lookup_2d_netcdf
   ! with that density class and filling the appropriate slice of the data array. Also
   ! checks for any errors, aborting if there were any.
   !
   ! Note that the lookup_2d_netcdf routine is designed to work with a single value of
   ! each of the indices. However, we want to fill parameter values for ALL density
   ! classes. This is why we loop over density class in this routine.
   !
   ! Note: inherits a number of variables from the parent routine

      use mkindexmapMod, only : lookup_2d_netcdf

      implicit none
      character(len=*), intent(in) :: varname       ! name of lookup table
      real(r8)        , intent(in) :: fill_val      ! value to put where we have no data in output variables
      logical         , intent(in) :: check_invalid ! should we check whether there are any invalid data in the output?
      real(r8)        , intent(out):: data(:,:)     ! output from lookup_2d_netcdf
      integer         , intent(in) :: n_extra_dims  ! number of extra dimensions in the lookup table

      ! slice to use if lookup table variable has more than 2 dimensions:
      type(dim_slice_type), intent(in), optional :: extra_dims(:)

      ! Local variables:

      integer :: k,n   ! indices
      integer :: ierr  ! error return code


      do k = 1, numurbl
         ! In the following, note that unity_dens_o(:,k) has been constructed so that
         ! unity_dens_o(:,k)==k everywhere. Thus, we fill data(:,k) with the parameter
         ! values corresponding to density class k.
         ! Also note: We use invalid_okay=.true. because we fill all density classes,
         ! some of which may have invalid entries. Because doing so disables some error
         ! checking, we do our own error checking after the call.
         call lookup_2d_netcdf(ncidi, varname, .true., &
                               'density_class', 'region', n_extra_dims, &
                               unity_dens_o(:,k), region_o, fill_val, data(:,k), ierr, &
                               extra_dims=extra_dims, nodata=index_nodata, &
                               invalid_okay=.true.)

         if (ierr /= 0) then
            write(6,*) modname//':'//subname//' ERROR in lookup_2d_netcdf for ', &
                 trim(varname), ' class', k, ': err=', ierr
            call abort()
         end if

         if (check_invalid) then
            ! Make sure we have valid parameter values wherever we have non-zero urban cover
            do n = 1, ns_o
               ! This check assumes that fill_val doesn't appear in any of the valid entries
               ! of the lookup table
               if (urbn_classes_gcell_o(n,k) > 0. .and. data(n,k) == fill_val) then
                  write(6,*) modname//':'//subname//' ERROR: fill value found in output where urban cover > 0'
                  write(6,*) 'var: ', trim(varname)
                  write(6,*) 'class: ', k
                  write(6,*) 'n: ', n
                  write(6,*) 'region: ', region_o(n)
                  write(6,*) 'urbn_classes_gcell_o(n,k): ', urbn_classes_gcell_o(n,k)
                  call abort()
               end if
            end do
         end if

      end do

   end subroutine lookup_and_check_err

end subroutine mkurbanpar
!------------------------------------------------------------------------------

end module mkurbanparMod
