module mkCH4inversionMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mkCH4inversionMod
!
! !DESCRIPTION:
! make inversion-derived parameters for CH4
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
  public mkCH4inversion            ! make inversion-derived parameters for CH4
!
!EOP
!===============================================================
contains
!===============================================================

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkCH4inversion
!
! !INTERFACE:
subroutine mkCH4inversion(ldomain, mapfname, datfname, ndiag, &
                          f0_o, p3_o, zwt0_o)
!
! !DESCRIPTION:
! make inversion-derived parameters for CH4
!
! !USES:
  use mkdomainMod, only : domain_type, domain_clean, domain_read
  use mkgridmapMod
  use mkncdio
  use mkdiagnosticsMod, only : output_diagnostics_continuous
  use mkchecksMod, only : min_bad, max_bad
!
! !ARGUMENTS:
  
  implicit none
  type(domain_type) , intent(in) :: ldomain
  character(len=*)  , intent(in) :: mapfname          ! input mapping file name
  character(len=*)  , intent(in) :: datfname          ! input data file name
  integer           , intent(in) :: ndiag             ! unit number for diag out
  real(r8)          , intent(out):: f0_o(:)           ! output grid: maximum gridcell fractional inundated area (unitless)
  real(r8)          , intent(out):: p3_o(:)           ! output grid: coefficient for qflx_surf_lag for finundated (s/mm)
  real(r8)          , intent(out):: zwt0_o(:)         ! output grid: decay factor for finundated (m)
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
  
  real(r8), parameter :: min_valid_f0    = 0._r8
  real(r8), parameter :: max_valid_f0    = 1._r8 + 1.0e-14_r8
  real(r8), parameter :: min_valid_p3    = 0._r8
  real(r8), parameter :: min_valid_zwt0  = 0._r8

  character(len=32) :: subname = 'mkCH4inversion'
!-----------------------------------------------------------------------

  write (6,*) 'Attempting to make inversion-derived CH4 parameters.....'
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

  write(6,*)'Open CH4 parameter file: ', trim(datfname)
  call check_ret(nf_open(datfname, 0, ncid), subname)

  allocate(data_i(tdomain%ns), stat=ier)
  if (ier/=0) call abort()

  ! -----------------------------------------------------------------
  ! Regrid f0
  ! -----------------------------------------------------------------

  call check_ret(nf_inq_varid (ncid, 'F0', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, data_i), subname)
  call gridmap_areaave(tgridmap, data_i, f0_o, nodata=0.01_r8)

  ! Check validity of output data
  if (min_bad(f0_o, min_valid_f0, 'f0') .or. &
      max_bad(f0_o, max_valid_f0, 'f0')) then
     stop
  end if

  call output_diagnostics_continuous(data_i, f0_o, tgridmap, "F0: max inundated area", "unitless", ndiag)

  ! -----------------------------------------------------------------
  ! Regrid p3
  ! -----------------------------------------------------------------

  call check_ret(nf_inq_varid (ncid, 'P3', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, data_i), subname)
  call gridmap_areaave(tgridmap, data_i, p3_o, nodata=10._r8)

  ! Check validity of output data
  if (min_bad(p3_o, min_valid_p3, 'p3')) then
     stop
  end if

  call output_diagnostics_continuous(data_i, p3_o, tgridmap, "P3: finundated coeff", "s/mm", ndiag)

  ! -----------------------------------------------------------------
  ! Regrid zwt0
  ! -----------------------------------------------------------------

  call check_ret(nf_inq_varid (ncid, 'ZWT0', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, data_i), subname)
  call gridmap_areaave(tgridmap, data_i, zwt0_o, nodata=0.01_r8)

  ! Check validity of output data
  if (min_bad(zwt0_o, min_valid_zwt0, 'zwt0')) then
     stop
  end if

  call output_diagnostics_continuous(data_i, zwt0_o, tgridmap, "ZWT0: finundated decay factor", "m", ndiag)

  ! -----------------------------------------------------------------
  ! Close files and deallocate dynamic memory
  ! -----------------------------------------------------------------

  call check_ret(nf_close(ncid), subname)
  call domain_clean(tdomain) 
  call gridmap_clean(tgridmap)
  deallocate (data_i)

  write (6,*) 'Successfully made inversion-derived CH4 parameters'
  write (6,*)
  call shr_sys_flush(6)

end subroutine mkCH4inversion


end module mkCH4inversionMod
