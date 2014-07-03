module mkVICparamsMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mkVICparamsMod
!
! !DESCRIPTION:
! make parameters for VIC
!
! !REVISION HISTORY:
! Author: Bill Sacks
!
!-----------------------------------------------------------------------
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_flush
  use mkdomainMod , only : domain_checksame

  implicit none

  private

! !PUBLIC MEMBER FUNCTIONS:
  public mkVICparams            ! make VIC parameters
!
!EOP
!===============================================================
contains
!===============================================================

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkVICparams
!
! !INTERFACE:
subroutine mkVICparams(ldomain, mapfname, datfname, ndiag, &
                       binfl_o, ws_o, dsmax_o, ds_o)
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
  real(r8)          , intent(out):: binfl_o(:)        ! output grid: VIC b parameter for the Variable Infiltration Capacity Curve (unitless)
  real(r8)          , intent(out):: ws_o(:)           ! output grid: VIC Ws parameter for the ARNO curve (unitless)
  real(r8)          , intent(out):: dsmax_o(:)        ! output grid: VIC Dsmax parameter for the ARNO curve (mm/day)
  real(r8)          , intent(out):: ds_o(:)           ! output grid: VIC Ds parameter for the ARNO curve (unitless)
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
  type(gridmap_type)    :: tgridmap
  type(domain_type)     :: tdomain            ! local domain
  real(r8), allocatable :: data_i(:)          ! data on input grid
  integer  :: ncid,varid                      ! input netCDF id's
  integer  :: ier                             ! error status
  
  real(r8), parameter :: min_valid_binfl = 0._r8
  real(r8), parameter :: min_valid_ws    = 0._r8
  real(r8), parameter :: min_valid_dsmax = 0._r8
  real(r8), parameter :: min_valid_ds    = 0._r8

  character(len=32) :: subname = 'mkVICparams'
!-----------------------------------------------------------------------

  write (6,*) 'Attempting to make VIC parameters.....'
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

  write(6,*)'Open VIC parameter file: ', trim(datfname)
  call check_ret(nf_open(datfname, 0, ncid), subname)

  allocate(data_i(tdomain%ns), stat=ier)
  if (ier/=0) call abort()

  ! -----------------------------------------------------------------
  ! Regrid binfl
  ! -----------------------------------------------------------------

  call check_ret(nf_inq_varid (ncid, 'binfl', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, data_i), subname)
  call gridmap_areaave(tgridmap, data_i, binfl_o, nodata=0.1_r8)

  ! Check validity of output data
  if (min_bad(binfl_o, min_valid_binfl, 'binfl')) then
     stop
  end if

  call output_diagnostics_continuous(data_i, binfl_o, tgridmap, "VIC b parameter", "unitless", ndiag)

  ! -----------------------------------------------------------------
  ! Regrid Ws
  ! -----------------------------------------------------------------

  call check_ret(nf_inq_varid (ncid, 'Ws', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, data_i), subname)
  call gridmap_areaave(tgridmap, data_i, ws_o, nodata=0.75_r8)

  ! Check validity of output data
  if (min_bad(ws_o, min_valid_ws, 'Ws')) then
     stop
  end if

  call output_diagnostics_continuous(data_i, ws_o, tgridmap, "VIC Ws parameter", "unitless", ndiag)

  ! -----------------------------------------------------------------
  ! Regrid Dsmax
  ! -----------------------------------------------------------------

  call check_ret(nf_inq_varid (ncid, 'Dsmax', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, data_i), subname)
  call gridmap_areaave(tgridmap, data_i, dsmax_o, nodata=10._r8)

  ! Check validity of output data
  if (min_bad(dsmax_o, min_valid_dsmax, 'Dsmax')) then
     stop
  end if

  call output_diagnostics_continuous(data_i, dsmax_o, tgridmap, "VIC Dsmax parameter", "mm/day", ndiag)

  ! -----------------------------------------------------------------
  ! Regrid Ds
  ! -----------------------------------------------------------------

  call check_ret(nf_inq_varid (ncid, 'Ds', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, data_i), subname)
  call gridmap_areaave(tgridmap, data_i, ds_o, nodata=0.1_r8)

  ! Check validity of output data
  if (min_bad(ds_o, min_valid_ds, 'Ds')) then
     stop
  end if

  call output_diagnostics_continuous(data_i, ds_o, tgridmap, "VIC Ds parameter", "unitless", ndiag)

  ! -----------------------------------------------------------------
  ! Close files and deallocate dynamic memory
  ! -----------------------------------------------------------------

  call check_ret(nf_close(ncid), subname)
  call domain_clean(tdomain) 
  call gridmap_clean(tgridmap)
  deallocate (data_i)

  write (6,*) 'Successfully made VIC parameters'
  write (6,*)
  call shr_sys_flush(6)

end subroutine mkVICparams


end module mkVICparamsMod
