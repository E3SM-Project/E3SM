module mkvocefMod
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mkvocMod
!
! !DESCRIPTION:
! Make VOC percentage emissions for surface dataset
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
!-----------------------------------------------------------------------
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_flush
  use mkdomainMod , only : domain_checksame

  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:
  public :: mkvocef  ! Get the percentage emissions for VOC for different
                     ! land cover types
!EOP

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkvocef
!
! !INTERFACE:
subroutine mkvocef(ldomain, mapfname, datfname, ndiag, &
                   ef_btr_o, ef_fet_o, ef_fdt_o, ef_shr_o, ef_grs_o, ef_crp_o)
!
! !DESCRIPTION:
! make volatile organic coumpunds (VOC) emission factors.
!
! !USES:
  use mkdomainMod, only : domain_type, domain_clean, domain_read
  use mkgridmapMod
  use mkvarpar	
  use mkvarctl    
  use mkncdio
!
! !ARGUMENTS:
  implicit none
  type(domain_type) , intent(in) :: ldomain
  character(len=*)  , intent(in) :: mapfname    ! input mapping file name
  character(len=*)  , intent(in) :: datfname    ! input data file name
  integer           , intent(in) :: ndiag       ! unit number for diagnostic output
  real(r8)          , intent(out):: ef_btr_o(:) ! output grid: EFs for broadleaf trees
  real(r8)          , intent(out):: ef_fet_o(:) ! output grid: EFs for fineleaf evergreen
  real(r8)          , intent(out):: ef_fdt_o(:) ! output grid: EFs for fineleaf deciduous
  real(r8)          , intent(out):: ef_shr_o(:) ! output grid: EFs for shrubs
  real(r8)          , intent(out):: ef_grs_o(:) ! output grid: EFs for grasses
  real(r8)          , intent(out):: ef_crp_o(:) ! output grid: EFs for crops
!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! Author: Colette L. Heald
! 17 Jul 2007 F Vitt -- updated to pftintdat06_clm3_5_05 and corrected indexing of ef_*_i arrarys
!
!EOP
!
! !LOCAL VARIABLES:
  type(gridmap_type)    :: tgridmap
  type(domain_type)     :: tdomain          ! local domain
  real(r8), allocatable :: ef_btr_i(:)      ! input grid: EFs for broadleaf trees
  real(r8), allocatable :: ef_fet_i(:)      ! input grid: EFs for fineleaf evergreen
  real(r8), allocatable :: ef_fdt_i(:)      ! input grid: EFs for fineleaf deciduous
  real(r8), allocatable :: ef_shr_i(:)      ! input grid: EFs for shrubs
  real(r8), allocatable :: ef_grs_i(:)      ! input grid: EFs for grasses
  real(r8), allocatable :: ef_crp_i(:)      ! input grid: EFs for crops
  real(r8) :: sum_fldo                      ! global sum of dummy input fld
  real(r8) :: sum_fldi                      ! global sum of dummy input fld
  integer  :: k,n,no,ni,ns_o,ns_i           ! indices
  integer  :: ncid,dimid,varid              ! input netCDF id's
  integer  :: ier                           ! error status
  real(r8) :: relerr = 0.00001_r8           ! max error: sum overlap wts ne 1
  character(len=32) :: subname = 'mkvocef'
!-----------------------------------------------------------------------

  write (6,*) 'Attempting to make VOC emission factors .....'
  call shr_sys_flush(6)

  ns_o = ldomain%ns

  ! -----------------------------------------------------------------
  ! Read input Emission Factors
  ! -----------------------------------------------------------------

  ! Obtain input grid info, read local fields

  call domain_read(tdomain,datfname)
  ns_i = tdomain%ns
  allocate(ef_btr_i(ns_i), ef_fet_i(ns_i), ef_fdt_i(ns_i), &
           ef_shr_i(ns_i), ef_grs_i(ns_i), ef_crp_i(ns_i), &
           stat=ier)
  if (ier/=0) call abort()

  write (6,*) 'Open VOC file: ', trim(datfname)
  call check_ret(nf_open(datfname, 0, ncid), subname)
  call check_ret(nf_inq_varid (ncid, 'ef_btr', varid), subname)
  call check_ret(nf_get_var_double(ncid, varid, ef_btr_i), subname)
  call check_ret(nf_inq_varid (ncid, 'ef_fet', varid), subname)
  call check_ret(nf_get_var_double(ncid, varid, ef_fet_i), subname)
  call check_ret(nf_inq_varid (ncid, 'ef_fdt', varid), subname)
  call check_ret(nf_get_var_double(ncid, varid, ef_fdt_i), subname)
  call check_ret(nf_inq_varid (ncid, 'ef_shr', varid), subname)
  call check_ret(nf_get_var_double(ncid, varid, ef_shr_i), subname)
  call check_ret(nf_inq_varid (ncid, 'ef_grs', varid), subname)
  call check_ret(nf_get_var_double(ncid, varid, ef_grs_i), subname)
  call check_ret(nf_inq_varid (ncid, 'ef_crp', varid), subname)
  call check_ret(nf_get_var_double(ncid, varid, ef_crp_i), subname)
  call check_ret(nf_close(ncid), subname)

  ! Area-average percent cover on input grid to output grid 
  ! and correct according to land landmask
  ! Note that percent cover is in terms of total grid area.

  call gridmap_mapread(tgridmap, mapfname )

  ! Error checks for domain and map consistencies

  call domain_checksame( tdomain, ldomain, tgridmap )

  ! Do mapping from input to output grid

  call gridmap_areaave(tgridmap, ef_btr_i, ef_btr_o, nodata=0._r8)
  call gridmap_areaave(tgridmap, ef_fet_i, ef_fet_o, nodata=0._r8)
  call gridmap_areaave(tgridmap, ef_fdt_i, ef_fdt_o, nodata=0._r8)
  call gridmap_areaave(tgridmap, ef_shr_i, ef_shr_o, nodata=0._r8)
  call gridmap_areaave(tgridmap, ef_grs_i, ef_grs_o, nodata=0._r8)
  call gridmap_areaave(tgridmap, ef_crp_i, ef_crp_o, nodata=0._r8)

  ! Check for conservation

  do no = 1, ns_o
     if ( ef_btr_o(no) < 0._r8 ) then
        write (6,*) 'MKVOCEF error: EF btr = ',ef_btr_o(no), &
             ' is negative for no = ',no
        call abort()
     end if
     if ( ef_fet_o(no) < 0._r8 ) then
        write (6,*) 'MKVOCEF error: EF fet = ',ef_fet_o(no), &
             ' is negative for no = ',no
        call abort()
     end if
     if ( ef_fdt_o(no) < 0._r8 ) then
        write (6,*) 'MKVOCEF error: EF fdt = ',ef_fdt_o(no), &
             ' is negative for no = ',no
        call abort()
     end if
     if ( ef_shr_o(no) < 0._r8 ) then
        write (6,*) 'MKVOCEF error: EF shr = ',ef_shr_o(no), &
             ' is negative for no = ',no
        call abort()
     end if
     if ( ef_grs_o(no) < 0._r8 ) then
        write (6,*) 'MKVOCEF error: EF grs = ',ef_grs_o(no), &
             ' is negative for no = ',no
        call abort()
     end if
     if ( ef_crp_o(no) < 0._r8 ) then
        write (6,*) 'MKVOCEF error: EF crp = ',ef_crp_o(no), &
             ' is negative for no = ',no
        call abort()
     end if
  enddo

  ! -----------------------------------------------------------------
  ! Error check1
  ! Compare global sum fld_o to global sum fld_i.
  ! -----------------------------------------------------------------

  ! Global sum of output field -- must multiply by fraction of
  ! output grid that is land as determined by input grid

  sum_fldi = 0.0_r8
  do ni = 1,ns_i
     sum_fldi = sum_fldi + tgridmap%area_src(ni) * tgridmap%frac_src(ni)
  enddo

  sum_fldo = 0._r8
  do no = 1,ns_o
     sum_fldo = sum_fldo + tgridmap%area_dst(no) * tgridmap%frac_dst(no)
  end do

  if ( trim(mksrf_gridtype) == 'global') then
     if ( abs(sum_fldo/sum_fldi-1._r8) > relerr ) then
        write (6,*) 'MKVOCEF error: input field not conserved'
        write (6,'(a30,e20.10)') 'global sum output field = ',sum_fldo
        write (6,'(a30,e20.10)') 'global sum input  field = ',sum_fldi
        stop
     end if
  end if

  write (6,*) 'Successfully made VOC Emission Factors'
  write (6,*)
  call shr_sys_flush(6)

  ! Deallocate dynamic memory

  deallocate ( ef_btr_i, ef_fet_i, ef_fdt_i, &
               ef_shr_i, ef_grs_i, ef_crp_i )
  call domain_clean(tdomain)
  call gridmap_clean(tgridmap)

end subroutine mkvocef

!-----------------------------------------------------------------------

end module mkvocefMod
