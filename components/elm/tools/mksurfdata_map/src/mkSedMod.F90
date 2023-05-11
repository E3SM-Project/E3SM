module mkSedMod
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mkSedMod
!
! !DESCRIPTION:
! Make soil erosion data (gravel, slope percentile and parameters)
!
! !REVISION HISTORY:
! Author: Zeli Tan
!
!-----------------------------------------------------------------------
!!USES:
  use shr_kind_mod, only : r8 => shr_kind_r8, r4=>shr_kind_r4
  use shr_sys_mod , only : shr_sys_flush
  use mkdomainMod , only : domain_checksame
  implicit none

  SAVE
  private           ! By default make data private
!
! !PUBLIC MEMBER FUNCTIONS:
!
  public mksedAtt      ! Add attributes to output file

  public mkgrvl         ! Set soil gravel
  public mkslp10        ! Set slope percentile
  public mkEROparams    ! Set ELM-Erosion parameters
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
! !IROUTINE: mkgrvl
!
! !INTERFACE:
subroutine mkgrvl(ldomain, mapfname, datfname, ndiag, grvl_o)
!
! !DESCRIPTION:
! make gravel dataset
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
  type(domain_type), intent(in) :: ldomain
  character(len=*)  , intent(in) :: mapfname       ! input mapping file name
  character(len=*)  , intent(in) :: datfname       ! input data file name
  integer           , intent(in) :: ndiag          ! unit number for diag out
  real(r8)          , intent(out):: grvl_o(:,:)    ! output grid:
!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! 
! Author: Zeli Tan
!
!
! !LOCAL VARIABLES:
!EOP
  type(gridmap_type)    :: tgridmap
  type(domain_type)    :: tdomain         ! local domain
  real(r8), allocatable :: grvl_i(:,:)     ! input grid: total column gravel 
  real(r8) :: sum_fldi                     ! global sum of dummy input fld
  real(r8) :: sum_fldo                     ! global sum of dummy output fld
  real(r8) :: gomlev_i                     ! input  grid: global gravel on lev
  real(r8) :: garea_i                      ! input  grid: global area
  real(r8) :: gomlev_o                     ! output grid: global gravel on lev
  real(r8) :: garea_o                      ! output grid: global area
  integer  :: k,n,m,ni,no,ns_i             ! indices
  integer  :: lev                          ! level index
  integer  :: nlay                         ! number of soil layers
  integer  :: ncid,dimid,varid             ! input netCDF id's
  integer  :: ier                          ! error status
  real(r8) :: relerr = 0.00001             ! max error: sum overlap wts ne 1
  character(len=32) :: subname = 'mkgrvl'
!-----------------------------------------------------------------------

  write (6,*) 'Attempting to make gravel dataset .....'
  call shr_sys_flush(6)

  ! -----------------------------------------------------------------
  ! Read input file
  ! -----------------------------------------------------------------

  ! Obtain input grid info, read local fields

  call domain_read(tdomain,datfname)
  ns_i = tdomain%ns

  write (6,*) 'Open gravel file: ', trim(datfname)
  call check_ret(nf_open(datfname, 0, ncid), subname)

  call check_ret(nf_inq_dimid  (ncid, 'number_of_layers', dimid), subname)
  call check_ret(nf_inq_dimlen (ncid, dimid, nlay), subname)

  allocate(grvl_i(ns_i,nlay),stat=ier)
  if (ier/=0) call abort()
  if (nlay /= nlevsoi) then
     write(6,*)'nlay, nlevsoi= ',nlay,nlevsoi,' do not match'
     stop
  end if

  call check_ret(nf_inq_varid (ncid, 'PCT_GRVL', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, grvl_i), subname)

  call check_ret(nf_close(ncid), subname)

  ! Area-average percent cover on input grid to output grid 
  ! and correct according to land landmask
  ! Note that percent cover is in terms of total grid area.

  call gridmap_mapread(tgridmap, mapfname )

  call domain_checksame( tdomain, ldomain, tgridmap )

  do lev = 1,nlay
     call gridmap_areaave(tgridmap, grvl_i(:,lev), grvl_o(:,lev), nodata=0._r8)
  end do

  do lev = 1,nlevsoi

     ! Check for conservation

     do no = 1,ldomain%ns
        if ((grvl_o(no,lev)) > 100.000001_r8) then
           write (6,*) 'MKGRVL error: grvl = ',grvl_o(no,lev), &
                ' greater than 100.000001 for column, row = ',no
           call shr_sys_flush(6)
           stop
        end if
     enddo

     call shr_sys_flush(ndiag)

     write (6,*) 'Successfully made gravel, level = ', lev
     call shr_sys_flush(6)

  end do   ! lev

  ! Deallocate dynamic memory

  call domain_clean(tdomain)
  call gridmap_clean(tgridmap)
  deallocate (grvl_i)

  write (6,*) 'Successfully made gravel'
  call shr_sys_flush(6)
  write(6,*)

end subroutine mkgrvl

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkslp10
!
! !INTERFACE:
subroutine mkslp10(ldomain, mapfname, datfname, ndiag, slp10_o)
!
! !DESCRIPTION:
! make slope percentile dataset
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
  type(domain_type), intent(in) :: ldomain
  character(len=*)  , intent(in) :: mapfname       ! input mapping file name
  character(len=*)  , intent(in) :: datfname       ! input data file name
  integer           , intent(in) :: ndiag          ! unit number for diag out
  real(r8)          , intent(out):: slp10_o(:,:)   ! output grid:
!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! 
! Author: Zeli Tan
!
!
! !LOCAL VARIABLES:
!EOP
  type(gridmap_type)    :: tgridmap
  type(domain_type)    :: tdomain         ! local domain
  real(r8), allocatable :: slp10_i(:,:)    ! input grid: total column slope
  real(r8) :: sum_fldi                     ! global sum of dummy input fld
  real(r8) :: sum_fldo                     ! global sum of dummy output fld
  real(r8) :: gomlev_i                     ! input  grid: global slope on lev
  real(r8) :: garea_i                      ! input  grid: global area
  real(r8) :: gomlev_o                     ! output grid: global slope on lev
  real(r8) :: garea_o                      ! output grid: global area
  integer  :: k,n,m,ni,no,ns_i             ! indices
  integer  :: lev                          ! level index
  integer  :: nlay                         ! number of slope percentile layers
  integer  :: ncid,dimid,varid             ! input netCDF id's
  integer  :: ier                          ! error status
  real(r8) :: relerr = 0.00001             ! max error: sum overlap wts ne 1
  character(len=32) :: subname = 'mkslp10'
!-----------------------------------------------------------------------

  write (6,*) 'Attempting to make slope percentile dataset .....'
  call shr_sys_flush(6)

  ! -----------------------------------------------------------------
  ! Read input file
  ! -----------------------------------------------------------------

  ! Obtain input grid info, read local fields

  call domain_read(tdomain,datfname)
  ns_i = tdomain%ns

  write (6,*) 'Open slope percentile file: ', trim(datfname)
  call check_ret(nf_open(datfname, 0, ncid), subname)

  call check_ret(nf_inq_dimid  (ncid, 'level', dimid), subname)
  call check_ret(nf_inq_dimlen (ncid, dimid, nlay), subname)

  allocate(slp10_i(ns_i,nlay),stat=ier)
  if (ier/=0) call abort()
  if (nlay /= nlevslp) then
     write(6,*)'nlay, nlevslp= ',nlay,nlevslp,' do not match'
     stop
  end if

  call check_ret(nf_inq_varid (ncid, 'SLP_P10', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, slp10_i), subname)

  call check_ret(nf_close(ncid), subname)

  ! Area-average percent cover on input grid to output grid 
  ! and correct according to land landmask
  ! Note that percent cover is in terms of total grid area.

  call gridmap_mapread(tgridmap, mapfname )

  call domain_checksame( tdomain, ldomain, tgridmap )

  do lev = 1,nlay
     call gridmap_areaave(tgridmap, slp10_i(:,lev), slp10_o(:,lev), nodata=0._r8)
  end do

  do lev = 1,nlevslp

     ! Check for conservation

     do no = 1,ldomain%ns
        if ((slp10_o(no,lev)) > 100.000001_r8) then
           write (6,*) 'MKGRVL error: slp10 = ',slp10_o(no,lev), &
                ' greater than 100.000001 for column, row = ',no
           call shr_sys_flush(6)
           stop
        end if
     enddo

     call shr_sys_flush(ndiag)

     write (6,*) 'Successfully made slope percentile, level = ', lev
     call shr_sys_flush(6)

  end do   ! lev

  ! Deallocate dynamic memory

  call domain_clean(tdomain)
  call gridmap_clean(tgridmap)
  deallocate (slp10_i)

  write (6,*) 'Successfully made slope percentile'
  call shr_sys_flush(6)
  write(6,*)

end subroutine mkslp10 

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkEROparams
!
! !INTERFACE:
subroutine mkEROparams(ldomain, mapfname, datfname, ndiag, &
                       ero_c1_o, ero_c2_o, ero_c3_o, tillage_o, &
                       litho_o)
!
! !DESCRIPTION:
! make VIC parameters
!
! !USES:
  use mkdomainMod, only : domain_type, domain_clean, domain_read
  use mkgridmapMod
  use mkncdio
  use mkdiagnosticsMod, only : output_diagnostics_continuous
  use mkchecksMod, only : min_bad
!
! !ARGUMENTS:

  implicit none
  type(domain_type) , intent(in) :: ldomain
  character(len=*)  , intent(in) :: mapfname          ! input mapping file name
  character(len=*)  , intent(in) :: datfname          ! input data file name
  integer           , intent(in) :: ndiag             ! unit number for diag out
  real(r8)          , intent(out):: ero_c1_o(:)       ! output grid: ELM-Erosion a scalar parameter for rainfall-driven hillslope erosion 
  real(r8)          , intent(out):: ero_c2_o(:)       ! output grid: ELM-Erosion a scalar parameter for runoff-driven hillslope erosion
  real(r8)          , intent(out):: ero_c3_o(:)       ! output grid: ELM-Erosion a scalar parameter for transport capacity of hillslope overland flow
  real(r8)          , intent(out):: tillage_o(:)      ! output grid: ELM-Erosion conserved tillage fraction
  real(r8)          , intent(out):: litho_o(:)        ! output grid: ELM-Erosion lithology erodibility index 
!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! Author: Zeli Tan
!
!
! !LOCAL VARIABLES:
!EOP
  type(gridmap_type)    :: tgridmap
  type(domain_type)     :: tdomain            ! local domain
  real(r8), allocatable :: data_i(:)          ! data on input grid
  integer  :: ncid,varid                      ! input netCDF id's
  integer  :: ier                             ! error status

  character(len=32) :: subname = 'mkEROparams'
!-----------------------------------------------------------------------

  write (6,*) 'Attempting to make ELM-Erosion parameters.....'
  call shr_sys_flush(6)

  ! -----------------------------------------------------------------
  ! Read domain and mapping information, check for consistency
  ! -----------------------------------------------------------------

  call domain_read(tdomain,datfname)

  call gridmap_mapread(tgridmap, mapfname )
  call gridmap_check( tgridmap, subname )

  call domain_checksame( tdomain, ldomain, tgridmap )

  ! -----------------------------------------------------------------
  ! Open input file, allocate memory for input data
  ! -----------------------------------------------------------------

  write(6,*)'Open ELM-Erosion parameter file: ', trim(datfname)
  call check_ret(nf_open(datfname, 0, ncid), subname)

  allocate(data_i(tdomain%ns), stat=ier)
  if (ier/=0) call abort()

  ! -----------------------------------------------------------------
  ! Regrid parEro_c1
  ! -----------------------------------------------------------------

  call check_ret(nf_inq_varid (ncid, 'parEro_c1', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, data_i), subname)
  call gridmap_areaave(tgridmap, data_i, ero_c1_o, nodata=0.0_r8)

  ! Check validity of output data
  if (min_bad(ero_c1_o, 0.0_r8, 'parEro_c1')) then
     stop
  end if

  call output_diagnostics_continuous(data_i, ero_c1_o, tgridmap, "ELM-Erosion c1 parameter", "unitless", ndiag)

  ! -----------------------------------------------------------------
  ! Regrid parEro_c2
  ! -----------------------------------------------------------------

  call check_ret(nf_inq_varid (ncid, 'parEro_c2', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, data_i), subname)
  call gridmap_areaave(tgridmap, data_i, ero_c2_o, nodata=0.0_r8)

  ! Check validity of output data
  if (min_bad(ero_c2_o, 0.0_r8, 'parEro_c2')) then
     stop
  end if

  call output_diagnostics_continuous(data_i, ero_c2_o, tgridmap, "ELM-Erosion c2 parameter", "unitless", ndiag)

  ! -----------------------------------------------------------------
  ! Regrid parEro_c3
  ! -----------------------------------------------------------------

  call check_ret(nf_inq_varid (ncid, 'parEro_c3', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, data_i), subname)
  call gridmap_areaave(tgridmap, data_i, ero_c3_o, nodata=0.0_r8)

  ! Check validity of output data
  if (min_bad(ero_c3_o, 0.0_r8, 'parEro_c3')) then
     stop
  end if

  call output_diagnostics_continuous(data_i, ero_c3_o, tgridmap, "ELM-Erosion c3 parameter", "unitless", ndiag)

  ! -----------------------------------------------------------------
  ! Regrid Tillage
  ! -----------------------------------------------------------------

  call check_ret(nf_inq_varid (ncid, 'Tillage', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, data_i), subname)
  call gridmap_areaave(tgridmap, data_i, tillage_o, nodata=0.0_r8)

  ! Check validity of output data
  if (min_bad(tillage_o, 0.0_r8, 'Tillage')) then
     stop
  end if

  call output_diagnostics_continuous(data_i, tillage_o, tgridmap, "ELM-Erosion CA fraction", "fraction", ndiag)

  ! -----------------------------------------------------------------
  ! Regrid Litho
  ! -----------------------------------------------------------------

  call check_ret(nf_inq_varid (ncid, 'Litho', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, data_i), subname)
  call gridmap_areaave(tgridmap, data_i, litho_o, nodata=0.0_r8)

  ! Check validity of output data
  if (min_bad(litho_o, 0.0_r8, 'Litho')) then
     stop
  end if

  call output_diagnostics_continuous(data_i, litho_o, tgridmap, "ELM-Erosion lithology erodibility", "unitless", ndiag)

  ! -----------------------------------------------------------------
  ! Close files and deallocate dynamic memory
  ! -----------------------------------------------------------------

  call check_ret(nf_close(ncid), subname)
  call domain_clean(tdomain)
  call gridmap_clean(tgridmap)
  deallocate (data_i)

  write (6,*) 'Successfully made ELM-Erosion parameters'
  write (6,*)
  call shr_sys_flush(6)

end subroutine mkEROparams

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mksedAtt
!
! !INTERFACE:
subroutine mksedAtt( ncid, dynlanduse, xtype )
!
! !DESCRIPTION:
! add atttributes to output file regarding the erosion module
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
! Original Author: Zeli Tan
!
!
! !LOCAL VARIABLES:
!EOP
  integer :: dimid                ! temporary
  character(len=256) :: str       ! global attribute string
  character(len=32) :: subname = 'mksedAtt'
!-----------------------------------------------------------------------

  if (.not. dynlanduse) then

     ! Define dimensions unique to erosion

     call check_ret(nf_def_dim (ncid, 'nlevslp',  &
                                       nlevslp    , dimid), subname)

     ! Add global attributes to file

     str = get_filename(mksrf_fgrvl)
     call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
          'soil_gravel_raw_data_file_name', len_trim(str), trim(str)), subname)

     str = get_filename(mksrf_fslp10)
     call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
          'slope_percentile_raw_data_file_name', len_trim(str), trim(str)), subname) 

     str = get_filename(mksrf_fero)
     call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
          'erosion_raw_data_file_name', len_trim(str), trim(str)), subname)
     
     ! Define variables

     if (outnc_1d) then
        call ncd_defvar(ncid=ncid, varname='PCT_GRVL', xtype=xtype, &
             dim1name='gridcell', dim2name='nlevsoi', &
             long_name='percent gravel', units='unitless')
     else
        call ncd_defvar(ncid=ncid, varname='PCT_GRVL', xtype=xtype, &
             dim1name='lsmlon', dim2name='lsmlat', dim3name='nlevsoi', &
             long_name='percent gravel', units='unitless')
     end if

     if (outnc_1d) then
        call ncd_defvar(ncid=ncid, varname='SLP_P10', xtype=xtype, &
             dim1name='gridcell', dim2name='nlevslp', &
             long_name='Slope at quantiles (minimum and 10 to 100 percentile)', &
             units='km km^-1')
     else
        call ncd_defvar(ncid=ncid, varname='SLP_P10', xtype=xtype, &
             dim1name='lsmlon', dim2name='lsmlat', dim3name='nlevslp', &
             long_name='Slope at quantiles (minimum and 10 to 100 percentile)', &
             units='km km^-1')
     end if

     if (outnc_1d) then
        call ncd_defvar(ncid=ncid, varname='parEro_c1', xtype=xtype, &
             dim1name='gridcell', & 
             long_name='a scalar parameter for rainfall-driven hillslope erosion', &
             units='unitless')
     else
        call ncd_defvar(ncid=ncid, varname='parEro_c1', xtype=xtype, &
             dim1name='lsmlon', dim2name='lsmlat', & 
             long_name='a scalar parameter for rainfall-driven hillslope erosion', &
             units='unitless')
     end if

     if (outnc_1d) then
        call ncd_defvar(ncid=ncid, varname='parEro_c2', xtype=xtype, &
             dim1name='gridcell', &
             long_name='a scalar parameter for runoff-driven hillslope erosion', &
             units='unitless')
     else
        call ncd_defvar(ncid=ncid, varname='parEro_c2', xtype=xtype, &
             dim1name='lsmlon', dim2name='lsmlat', &
             long_name='a scalar parameter for runoff-driven hillslope erosion', &
             units='unitless')
     end if

     if (outnc_1d) then
        call ncd_defvar(ncid=ncid, varname='parEro_c3', xtype=xtype, &
             dim1name='gridcell', &
             long_name='a scalar parameter for transport capacity of hillslope overland flow', &
             units='unitless')
     else    
        call ncd_defvar(ncid=ncid, varname='parEro_c3', xtype=xtype, &
             dim1name='lsmlon', dim2name='lsmlat', &
             long_name='a scalar parameter for transport capacity of hillslope overland flow', &
             units='unitless')
     end if

     if (outnc_1d) then
        call ncd_defvar(ncid=ncid, varname='Tillage', xtype=xtype, &
             dim1name='gridcell', &
             long_name='conserved tillage fraction', &
             units='fraction')
     else
        call ncd_defvar(ncid=ncid, varname='Tillage', xtype=xtype, &
             dim1name='lsmlon', dim2name='lsmlat', &
             long_name='conserved tillage fraction', &
             units='fraction')
     end if

     if (outnc_1d) then
        call ncd_defvar(ncid=ncid, varname='Litho', xtype=xtype, &
             dim1name='gridcell', &
             long_name='lithology erodibility index', &
             units='unitless')
     else
        call ncd_defvar(ncid=ncid, varname='Litho', xtype=xtype, &
             dim1name='lsmlon', dim2name='lsmlat', &
             long_name='lithology erodibility index', &
             units='unitless')
     end if

  end if

end subroutine mksedAtt

!-----------------------------------------------------------------------

end module mkSedMod
