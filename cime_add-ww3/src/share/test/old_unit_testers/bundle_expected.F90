!
! Module to test that a bundles is the same (or close to) the expected bundle
!
module bundle_expected
use shr_kind_mod
use shr_date_mod
use dshr_kind
use dshr_bundle
use dshr_domain

implicit none

private

public bundle_is_expected
public bundle_closeto_expected
public bundle_fill_cosz

private bundle_metadata_is_expected

contains

logical function bundle_is_expected( bun, expected_bun )
  use dshr_bundle, only : dshr_bundle_assignPtr
  implicit none

  type(dshr_bundle_bundleType), intent(IN) :: bun           ! bundle to test
  type(dshr_bundle_bundleType), intent(IN) :: expected_bun  ! expected bundle

  logical :: status
  real(r8), pointer :: data(:,:,:)
  real(r8), pointer :: exp_data(:,:,:)

  status = bundle_metadata_is_expected( bun, expected_bun )
  ! Compare field data
  if ( status )then
     call dshr_bundle_assignPtr( bun, data )
     call dshr_bundle_assignPtr( expected_bun, exp_data )
     if ( any(data /= exp_data) ) status = .false.
  end if
  bundle_is_expected = status

end function bundle_is_expected

logical function bundle_closeto_expected( bun, expected_bun, eps, crit_type )
  use dshr_bundle, only : dshr_bundle_assignPtr, dshr_bundle_getDims
  use shr_sys_mod
  implicit none
  type(dshr_bundle_bundleType), intent(IN) :: bun           ! bundle to test
  type(dshr_bundle_bundleType), intent(IN) :: expected_bun  ! expected bundle
  real(R8),                     intent(IN) :: eps           ! how close needs to be
  character(len=*), optional,   intent(IN) :: crit_type     ! criterea to determine if close enough

  type(dshr_domain_domainType),pointer :: domain
  real(r8), pointer :: data(:,:,:)
  real(r8), pointer :: exp_data(:,:,:), mask(:,:)
  real(r8) :: diff, max_diff, meansq, rms_diff
  integer :: i, j, f, ni, nj, nf, ndata
  character(len=90) :: lcrittype
  character(len=*), parameter :: subname = "bundle_closeto_expected"

  lcrittype = "abs_diff"
  if ( present(crit_type) ) lcrittype = crit_type
  if ( trim(lcrittype) /= "abs_diff" .and. trim(lcrittype) /= "rms_diff" .and. trim(lcrittype) /= "rel_diff" ) then
      call shr_sys_abort( subname//":: ERROR: bad criterea type input" )
  end if
  bundle_closeto_expected = bundle_metadata_is_expected( bun, expected_bun )
  ! Compare field data
  if ( bundle_closeto_expected )then
         call dshr_bundle_domainPtr( bun, domain )
         call dshr_bundle_getDims (bun,ni,nj,nf)
         allocate( mask(ni,nj) )
         call dshr_domain_getData( domain, mask, "mask" )
         call dshr_bundle_assignPtr( bun, data )
         call dshr_bundle_assignPtr( expected_bun, exp_data )
         max_diff = 0.0_r8
         meansq   = 0.0_r8
         ndata    = 0
outloop: do f = 1, nf
         do j = 1, nj
         do i = 1, ni
           if ( mask(i,j)  > 0.0_r8 )then
              diff = abs(data(i,j,f) - exp_data(i,j,f))
              ndata = ndata + 1
              meansq = meansq + diff**2
              if ( trim(lcrittype) == "rel_diff" .and. diff > 0.0_r8 ) diff = diff / max( abs(data(i,j,f)), abs(exp_data(i,j,f)) )
              if ( diff > max_diff ) max_diff = diff
              if ( diff > eps .and. .not. trim(lcrittype) == "rms_diff" )then
                 bundle_closeto_expected = .false.
              end if
           end if
         end do
         end do
         end do outloop
         deallocate( mask )
         rms_diff = sqrt(meansq/ndata)
         if ( rms_diff > eps .and. trim(lcrittype) == "rms_diff" ) bundle_closeto_expected = .false.
         write(*,*) "bundle_closeto_expected:  max_diff = ", max_diff, " RMS diff = ", rms_diff
  end if

end function bundle_closeto_expected

logical function bundle_metadata_is_expected( bun, expected_bun )
   use dshr_bundle, only : dshr_bundle_domainPtr, dshr_bundle_getDims, dshr_bundle_getFieldList, &
                           dshr_bundle_getDate
   use dshr_domain, only : dshr_domain_compare
   implicit none

   type(dshr_bundle_bundleType), intent(IN) :: bun           ! bundle to test
   type(dshr_bundle_bundleType), intent(IN) :: expected_bun  ! expected bundle

   type(dshr_domain_domainType),pointer :: domain
   type(dshr_domain_domainType),pointer :: exp_domain
   logical :: status
   integer :: ni, nj, nf, exp_ni, exp_nj, exp_nf
   integer :: date, sec, exp_date, exp_sec
   character(SHR_KIND_CX) :: fldlist, exp_fldlist


   call dshr_bundle_domainPtr( bun, domain )
   call dshr_bundle_domainPtr( expected_bun, exp_domain )

   status = dshr_domain_compare(    domain, exp_domain, method=dshr_domain_compareMaskIdent, eps=0.0_r8 )
   if ( status )then
      status = dshr_domain_compare( domain, exp_domain, method=dshr_domain_compareXYabs,     eps=0.0_r8 )
   end if
   if ( status )then
      call dshr_bundle_getDims( bun, ni, nj, nf )
      call dshr_bundle_getDims( bun, exp_ni, exp_nj, exp_nf )
      if ( ni /= exp_ni .or. nj /= exp_nj .or. nf /= exp_nf ) status = .false.
   end if
   if ( status )then
      call dshr_bundle_getFieldList( bun, fldlist )
      call dshr_bundle_getFieldList( expected_bun, exp_fldlist )
      if ( trim(fldlist) /= trim(exp_fldlist) ) status = .false.
   end if
   if ( status )then
      call dshr_bundle_getDate (bun,date,sec)
      call dshr_bundle_getDate (expected_bun,exp_date,exp_sec)
      if ( date /= exp_date .or. sec /= exp_sec ) status = .false.
   end if

   bundle_metadata_is_expected = status

end function bundle_metadata_is_expected

subroutine bundle_fill_cosz( scale, orb_eccen, orb_mvelpp, orb_lambm0, orb_obliqr, sdate_ub, domain, bun, kfld )
! Fill a bundle with data scaled by the average cosine of the solar zenith angle
   use shr_string_mod
   use shr_const_mod
   use shr_orb_mod
   use dshr_domain
   use shr_sys_mod
   implicit none
   real(r8), intent(IN) :: scale
   real(r8), intent(IN) :: orb_eccen, orb_mvelpp, orb_lambm0, orb_obliqr
   type(shr_date), intent(IN) :: sdate_ub             ! Upper bound of date for
   type(dshr_domain_domainType), pointer :: domain
   type(dshr_bundle_bundleType), intent(INOUT) :: bun ! bundle to fill
   integer, intent(in) :: kfld    ! Which field number to fill

   character(len=*), parameter :: subname = "bundle_fill_cosz"
   real(r8), pointer :: data(:,:,:), lat(:,:), lon(:,:), sumcosz(:,:)
   real(r8) :: cosz, calday, declin, eccf, calday_end
   integer :: i, j, f, ni, nj, nf, rc, t, ntimes, date_lb, sec_lb
   integer, parameter :: dtime = 18
   type(shr_date) :: sdate

   call dshr_domain_getDims(domain,ni,nj)
   allocate( lat(ni,nj) )
   allocate( lon(ni,nj) )
   allocate( sumcosz(ni,nj) )
   call dshr_domain_getData( domain, lat, "lat" )
   call dshr_domain_getData( domain, lon, "lon" )
   lat = lat * SHR_CONST_PI / 180._r8
   lon = lon * SHR_CONST_PI / 180._r8

   call dshr_bundle_assignPtr( bun, data )
   call dshr_bundle_getDate( bun, cdate=date_lb, sec=sec_lb )
   sdate = shr_date_initCDate( date_lb, 3600*24/dtime, sec_lb )
   calday_end = shr_date_getJulian( sdate_ub )
   sumcosz(:,:) = 0.0_r8
   calday = 0.0_r8
   ntimes = 0
   calday = shr_date_getJulian( sdate )
   nf = size( data, 3 )
   if ( kfld <= 0 .or. kfld > nf ) call shr_sys_abort( 'input kfld is out of bounds' )
   do while( sdate < sdate_ub .or. sdate == sdate_ub )
      ntimes = ntimes + 1
      call shr_orb_decl(calday ,orb_eccen ,orb_mvelpp ,orb_lambm0 ,orb_obliqr ,declin,eccf)
      do j = 1, nj
      do i = 1, ni
         cosz = shr_orb_cosz(calday,lat(i,j),lon(i,j),declin)
         if ( cosz < 0.01_r8 )  cosz = 0.01_r8
         if ( cosz < 0.001_r8 ) cosz = 0.001_r8
         sumcosz(i,j) = cosz + sumcosz(i,j)
      end do
      end do
      call shr_date_adv1step( sdate )
      calday = shr_date_getJulian( sdate )
   end do
   data(:,:,kfld) = sumcosz(:,:)*scale/real(ntimes,r8)

   nullify( data )
   nullify( domain )
   deallocate( lat )
   deallocate( lon )
   deallocate( sumcosz )

end subroutine bundle_fill_cosz

end module bundle_expected
