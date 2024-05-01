module mktopradMod
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mktopradMod
!
! !DESCRIPTION:
! Make topography data for TOP solar radiation parameterization
!
! !REVISION HISTORY:
! Author: Dalei Hao
!
!-----------------------------------------------------------------------
!!USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_flush
  implicit none

  SAVE
  private           ! By default make data private
!
! !PUBLIC MEMBER FUNCTIONS:
!
  public mktopradAtt      ! Add attributes to output file

  public mktoprad         ! Set topography
!
! !PUBLIC DATA MEMBERS:
!
!
! !PRIVATE DATA MEMBERS:
!
! !PRIVATE MEMBER FUNCTIONS:

!EOP
!===============================================================
contains
!===============================================================

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mktoprad
!
! !INTERFACE:
subroutine mktoprad(ldomain, mapfname, datfname, varname, ndiag, top_o, nodata)
!
! !DESCRIPTION:
! Make topography data for TOP solar radiation parameterization
!
! !USES:
  use mkdomainMod  , only : domain_type, domain_clean, domain_read, domain_checksame
  use mkgridmapMod
  use mkvarpar
  use mkvarctl
  use mkncdio
!
! !ARGUMENTS:
  implicit none
  type(domain_type), intent(in) :: ldomain
  character(len=*)  , intent(in) :: mapfname  ! input mapping file name
  character(len=*)  , intent(in) :: datfname  ! input data file name
  integer           , intent(in) :: ndiag     ! unit number for diag out
  character(len=*)  , intent(in) :: varname   ! topo variable name
  real(r8)          , intent(out):: top_o(:)  ! output topography data
  real(r8)          , intent(in) :: nodata    ! default value
!
!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! Author: Dalei Hao
!
!
! !LOCAL VARIABLES:
!EOP
  type(domain_type)     :: tdomain            ! local domain
  type(gridmap_type)    :: tgridmap           ! local gridmap

  real(r8), allocatable :: top_i(:)           ! input top variable
  real(r8), allocatable :: mask_i(:)          ! input grid: mask (0, 1)
  integer  :: ns_i,ns_o                       ! indices
  integer  :: k,l,n,m,ni                      ! indices
  integer  :: ncidi,dimid,varid               ! input netCDF id's
  integer  :: ier                             ! error status
  character(len=256) :: name                  ! name of attribute
  character(len=256) :: unit                  ! units of attribute
  character(len= 32) :: subname = 'mktop'
!-----------------------------------------------------------------------

  write (6,*) 'Attempting to make topography .....'
  call shr_sys_flush(6)

  ns_o = ldomain%ns

  ! -----------------------------------------------------------------
  ! Read input file
  ! -----------------------------------------------------------------

  ! Obtain input grid info, read local fields

  call domain_read(tdomain,datfname)

  ns_i = tdomain%ns
  allocate(top_i(ns_i), stat=ier)
  if (ier /= 0) then
     write(6,*)'mktoprad allocation error'; call abort()
  end if

  write (6,*) 'Open topography file: ', trim(datfname)
  call check_ret(nf_open(datfname, 0, ncidi), subname)
  call check_ret(nf_inq_varid (ncidi, trim(varname), varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, top_i), subname)
  call check_ret(nf_close(ncidi), subname)

  ! set mask as 0 when topo data is filled value: -9999
  allocate(mask_i(ns_i), stat=ier)
  if (ier /= 0) then
     write(6,*)'mktoprad allocation error'; call abort()
  end if
  
  mask_i(:) = 1._r8
  do ni = 1,ns_i
      if (top_i(ni) < -1000._r8) then
         mask_i(ni) = 0._r8
     end if
  enddo

  ! Read mapping file
  call gridmap_mapread(tgridmap, mapfname)

  ! Error checks for domain and map consistencies
  call domain_checksame( tdomain, ldomain, tgridmap )

  ! Determine top_o on output grid
  top_o(:) = nodata

  call gridmap_areaave(tgridmap, top_i, top_o, nodata=nodata, mask_src=mask_i)

  ! Deallocate dynamic memory

  call domain_clean(tdomain)
  call gridmap_clean(tgridmap)
  deallocate (top_i)
  deallocate (mask_i)

  write (6,*) 'Successfully made topography parameters'
  write (6,*)
  call shr_sys_flush(6)

end subroutine mktoprad

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mktopradAtt
!
! !INTERFACE:
subroutine mktopradAtt( ncid, dynlanduse, xtype )
!
! !DESCRIPTION:
! add atttributes to output file regarding the topography module
!
! !USES:
  use fileutils  , only : get_filename
  use mkncdio    , only : check_ret, ncd_defvar
  use mkvarpar   
  use mkvarctl   

! !ARGUMENTS:
  implicit none
  include 'netcdf.inc'
  integer, intent(in) :: ncid         ! NetCDF file ID to write out to
  logical, intent(in) :: dynlanduse   ! if dynamic land-use file
  integer, intent(in) :: xtype        ! external type to output real data as
!
! !CALLED FROM:
! subroutine mkfile in module mkfileMod
!
! !REVISION HISTORY:
! Original Author: Dalei Hao
!
!
! !LOCAL VARIABLES:
!EOP
  integer :: dimid                ! temporary
  character(len=256) :: str       ! global attribute string
  character(len=32) :: subname = 'mktopAtt'
!-----------------------------------------------------------------------

  if (.not. dynlanduse) then

     ! Add global attributes to file

     str = get_filename(mksrf_fgrvl)
     call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
          'top_raw_data_file_name', len_trim(str), trim(str)), subname)
     
     ! Define variables

     if (outnc_1d) then
        call ncd_defvar(ncid=ncid, varname='SINSL_COSAS', xtype=xtype, &
             dim1name='gridcell',&
             long_name='sin(slope) * cos(aspect)', units='unitless')
     else
        call ncd_defvar(ncid=ncid, varname='SINSL_COSAS', xtype=xtype, &
             dim1name='lsmlon', dim2name='lsmlat', &
             long_name='sin(slope) * cos(aspect)', units='unitless')
     end if

     if (outnc_1d) then
        call ncd_defvar(ncid=ncid, varname='SINSL_SINAS', xtype=xtype, &
             dim1name='gridcell',&
             long_name='sin(slope) * sin(aspect)', units='unitless')
     else
        call ncd_defvar(ncid=ncid, varname='SINSL_SINAS', xtype=xtype, &
             dim1name='lsmlon', dim2name='lsmlat', &
             long_name='sin(slope) * sin(aspect)', units='unitless')
     end if

    if (outnc_1d) then
        call ncd_defvar(ncid=ncid, varname='SKY_VIEW', xtype=xtype, &
            dim1name='gridcell',&
            long_name='sky view factor', units='unitless')
    else
        call ncd_defvar(ncid=ncid, varname='SKY_VIEW', xtype=xtype, &
            dim1name='lsmlon', dim2name='lsmlat', &
            long_name='sky view factor', units='unitless')
    end if

    if (outnc_1d) then
        call ncd_defvar(ncid=ncid, varname='TERRAIN_CONFIG', xtype=xtype, &
            dim1name='gridcell',&
            long_name='terrain configuration factor', units='unitless')
    else
        call ncd_defvar(ncid=ncid, varname='TERRAIN_CONFIG', xtype=xtype, &
            dim1name='lsmlon', dim2name='lsmlat', &
            long_name='terrain configuration factor', units='unitless')
    end if

  end if

end subroutine mktopradAtt

!-----------------------------------------------------------------------

end module mktopradMod
