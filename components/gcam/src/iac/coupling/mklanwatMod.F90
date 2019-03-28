module mklanwatMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mklanwatMod
!
! !DESCRIPTION:
! make %lake and %wetland from Cogley's one degree data
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!-----------------------------------------------------------------------
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_flush

  implicit none

  private

! !PUBLIC MEMBER FUNCTIONS:
  public mklanwat

!EOP
!===============================================================
contains
!===============================================================

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mklanwat
!
! !INTERFACE:
subroutine mklanwat(ldomain, mapfname, datfname, ndiag, zero_out, lake_o, swmp_o)
!
! !DESCRIPTION:
! make %lake and %wetland from Cogley's one degree data
!
! !USES:
  use mkfileutils , only : getfil
  use mkdomainMod, only : domain1_type, domain1_clean, domain1_read
  use mkgridmapMod
  use mkvarpar	
  use mkvarctl    
  use mkncdio
!
! !ARGUMENTS:
  
  implicit none
  type(domain1_type), intent(in) :: ldomain
  character(len=*)  , intent(in) :: mapfname  ! input mapping file name
  character(len=*)  , intent(in) :: datfname  ! input data file name
  integer           , intent(in) :: ndiag     ! unit number for diag out
  logical           , intent(in) :: zero_out  ! if should zero glacier out
  real(r8)          , intent(out):: lake_o(:) ! output grid: %lake
  real(r8)          , intent(out):: swmp_o(:) ! output grid: %wetland
!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
  type(gridmap_type)    :: tgridmap
  type(domain1_type)    :: tdomain            ! local domain
  real(r8), allocatable :: lake_i(:)          ! input grid: percent lake
  real(r8), allocatable :: swmp_i(:)          ! input grid: percent swamp
  real(r8) :: sum_fldi                        ! global sum of dummy input fld
  real(r8) :: sum_fldo                        ! global sum of dummy output fld
  real(r8) :: glake_i                         ! input  grid: global lake
  real(r8) :: gswmp_i                         ! input  grid: global swamp
  real(r8) :: garea_i                         ! input  grid: global area
  real(r8) :: glake_o                         ! output grid: global lake
  real(r8) :: gswmp_o                         ! output grid: global swamp
  real(r8) :: garea_o                         ! output grid: global area
  integer  :: ni,no,k,n,m,ns_i,ns_o           ! indices
  integer  :: ncid,dimid,varid                ! input netCDF id's
  integer  :: ier                             ! error status
  real(r8) :: relerr = 0.00001                ! max error: sum overlap wts ne 1
  character(len=256) locfn                    ! local dataset file name
  character(len=32) :: subname = 'mklanwat'
!-----------------------------------------------------------------------

  write (6,*) 'Attempting to make %lake and %wetland .....'
  call shr_sys_flush(6)

  ! -----------------------------------------------------------------
  ! Read input file
  ! -----------------------------------------------------------------

  ! Obtain input grid info, read local fields

  ns_o = ldomain%ns

  call getfil (datfname, locfn, 0)

  call domain1_read(tdomain,locfn)  !todo - put in routine to read in 1d domains 
  ns_i = tdomain%ns
  allocate(lake_i(ns_i), swmp_i(ns_i), stat=ier)
  if (ier/=0) call abort()

  call check_ret(nf_open(locfn, 0, ncid), subname)
  call check_ret(nf_inq_varid (ncid, 'PCT_LAKE', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, lake_i), subname)
  call check_ret(nf_inq_varid (ncid, 'PCT_WETLAND', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, swmp_i), subname)
  call check_ret(nf_close(ncid), subname)

  ! Area-average percent cover on input grid to output grid 
  ! and correct according to land landmask
  ! Note that percent cover is in terms of total grid area.

  call gridmap_mapread(tgridmap, mapfname )

  ! Error checks for domain and map consistencies

  if (tdomain%ns /= tgridmap%na) then
     write(6,*)'input domain size and gridmap source size are not the same size'
     write(6,*)' domain size = ',tdomain%ns
     write(6,*)' map src size= ',tgridmap%na
     stop
  end if
  do n = 1,tgridmap%ns
     ni = tgridmap%src_indx(n)
     if (tdomain%mask(ni) /= tgridmap%mask_src(ni)) then
        write(6,*)'input domain mask and gridmap mask are not the same at ni = ',ni
        write(6,*)' domain  mask= ',tdomain%mask(ni)
        write(6,*)' gridmap mask= ',tgridmap%mask_src(ni)
        stop
     end if
     if (tdomain%lonc(ni) /= tgridmap%xc_src(ni)) then
        write(6,*)'input domain lon and gridmap lon not the same at ni = ',ni
        write(6,*)' domain  lon= ',tdomain%lonc(ni)
        write(6,*)' gridmap lon= ',tgridmap%xc_src(ni)
        stop
     end if
     if (tdomain%latc(ni) /= tgridmap%yc_src(ni)) then
        write(6,*)'input domain lat and gridmap lat not the same at ni = ',ni
        write(6,*)' domain  lat= ',tdomain%latc(ni)
        write(6,*)' gridmap lat= ',tgridmap%yc_src(ni)
        stop
     end if
  end do

  ! Determine lake_o and swmp_o on output grid

  call gridmap_areaave(tgridmap, lake_i,lake_o)
  call gridmap_areaave(tgridmap, swmp_i,swmp_o)

  do no = 1,ns_o
     if (lake_o(no) < 1.) lake_o(no) = 0.
     if (swmp_o(no) < 1.) swmp_o(no) = 0.
     if (zero_out       ) lake_o(no) = 0.
     if (zero_out       ) swmp_o(no) = 0.
  enddo

  ! Check for conservation

  do no = 1,ns_o
     if ((lake_o(no) + swmp_o(no)) > 100.000001_r8) then
        write (6,*) 'MKLANWAT error: lake = ',lake_o(no), &
             ' and wetland = ',swmp_o(no), &
             ' sum are greater than 100 for no = ',no
        stop
     end if
  enddo

  ! -----------------------------------------------------------------
  ! Error check prep
  ! Global sum of output field -- must multiply by fraction of
  ! output grid that is land as determined by input grid
  ! -----------------------------------------------------------------

  sum_fldi = 0.0_r8
  do ni = 1,ns_i
    sum_fldi = sum_fldi + tgridmap%area_src(ni)*tgridmap%frac_src(ni)*re**2
  enddo

  sum_fldo = 0.
  do no = 1,ns_o
     sum_fldo = sum_fldo + tgridmap%area_dst(no)*tgridmap%frac_dst(no)*re**2
  end do

  ! -----------------------------------------------------------------
  ! Error check1
  ! Compare global sum fld_o to global sum fld_i.
  ! -----------------------------------------------------------------

  if ( .not. zero_out .and. (trim(mksrf_gridtype) == 'global') ) then
     if ( abs(sum_fldo/sum_fldi-1.) > relerr ) then
        write (6,*) 'MKLANWAT error: input field not conserved'
        write (6,'(a30,e20.10)') 'global sum output field = ',sum_fldo
        write (6,'(a30,e20.10)') 'global sum input  field = ',sum_fldi
        stop
     end if
  end if

  ! -----------------------------------------------------------------
  ! Error check2
  ! Compare global areas on input and output grids
  ! -----------------------------------------------------------------

  ! Input grid

  glake_i = 0.
  gswmp_i = 0.
  garea_i = 0.
  do ni = 1,ns_i
     garea_i = garea_i + tgridmap%area_src(ni)*re**2
     glake_i = glake_i + lake_i(ni)*tgridmap%area_src(ni)/100.*re**2
     gswmp_i = gswmp_i + swmp_i(ni)*tgridmap%area_src(ni)/100.*re**2
  end do

  ! Output grid

  glake_o = 0.
  gswmp_o = 0.
  garea_o = 0.
  do no = 1,ns_o
     garea_o = garea_o + tgridmap%area_dst(no)*re**2
     glake_o = glake_o + lake_o(no)*tgridmap%area_dst(no)/100.*re**2
     gswmp_o = gswmp_o + swmp_o(no)*tgridmap%area_dst(no)/100.*re**2
  end do

  ! Diagnostic output

  write (ndiag,*)
  write (ndiag,'(1x,70a1)') ('=',k=1,70)
  write (ndiag,*) 'Inland Water Output'
  write (ndiag,'(1x,70a1)') ('=',k=1,70)

  write (ndiag,*)
  write (ndiag,'(1x,70a1)') ('.',k=1,70)
  write (ndiag,2001)
2001 format (1x,'surface type   input grid area  output grid area'/ &
             1x,'                 10**6 km**2      10**6 km**2   ')
  write (ndiag,'(1x,70a1)') ('.',k=1,70)
  write (ndiag,*)
  write (ndiag,2002) glake_i*1.e-06,glake_o*1.e-06
  write (ndiag,2003) gswmp_i*1.e-06,gswmp_o*1.e-06
  write (ndiag,2004) garea_i*1.e-06,garea_o*1.e-06
2002 format (1x,'lakes       ',f14.3,f17.3)
2003 format (1x,'wetlands    ',f14.3,f17.3)
2004 format (1x,'all surface ',f14.3,f17.3)

  ! Deallocate dynamic memory

  call domain1_clean(tdomain) 
  call gridmap_clean(tgridmap)
  deallocate (lake_i,swmp_i)

  write (6,*) 'Successfully made %lake and %wetland'
  write (6,*)
  call shr_sys_flush(6)

end subroutine mklanwat

end module mklanwatMod
