module mkdiagnosticsMod
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mkdiagnostics
!
! !DESCRIPTION:
! Output diagnostics to log file
!
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  
  implicit none
  private
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: output_diagnostics_area               ! output diagnostics for field that is % of grid area
  public :: output_diagnostics_continuous         ! output diagnostics for a continuous (real-valued) field
  public :: output_diagnostics_continuous_outonly ! output diagnostics for a continuous (real-valued) field, just on the output grid
!
!
! !REVISION HISTORY:
! Author: Bill Sacks
!
!EOP
!------------------------------------------------------------------------------
contains

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: output_diagnostics_area
!
! !INTERFACE:
subroutine output_diagnostics_area(data_i, data_o, gridmap, name, percent, ndiag)
!
! !DESCRIPTION:
! Output diagnostics for a field that gives either fraction or percent of grid cell area
!
! !USES:
  use mkgridmapMod, only : gridmap_type
  use mkvarpar, only : re
!
! !ARGUMENTS:
  implicit none
  real(r8)          , intent(in) :: data_i(:)    ! data on input grid
  real(r8)          , intent(in) :: data_o(:)    ! data on output grid
  type(gridmap_type), intent(in) :: gridmap      ! mapping info
  character(len=*)  , intent(in) :: name         ! name of field
  logical           , intent(in) :: percent      ! is field specified as percent? (alternative is fraction)
  integer           , intent(in) :: ndiag        ! unit number for diagnostic output
!
! !REVISION HISTORY:
! Author: Bill Sacks
!
!
! !LOCAL VARIABLES:
!EOP
  real(r8) :: gdata_i         ! global sum of input data
  real(r8) :: gdata_o         ! global sum of output data
  real(r8) :: garea_i         ! global sum of input area
  real(r8) :: garea_o         ! global sum of output area
  integer  :: ns_i, ns_o      ! sizes of input & output grids
  integer  :: ni,no,k         ! indices

  character(len=*), parameter :: subname = "output_diagnostics_area"
!------------------------------------------------------------------------------

  ! Error check for array size consistencies

  ns_i = gridmap%na
  ns_o = gridmap%nb
  if (size(data_i) /= ns_i .or. &
      size(data_o) /= ns_o) then
     write(6,*) subname//' ERROR: array size inconsistencies for ', trim(name)
     write(6,*) 'size(data_i) = ', size(data_i)
     write(6,*) 'ns_i         = ', ns_i
     write(6,*) 'size(data_o) = ', size(data_o)
     write(6,*) 'ns_o         = ', ns_o
     stop
  end if

  ! Sums on input grid

  gdata_i = 0.
  garea_i = 0.
  do ni = 1,ns_i
     garea_i = garea_i + gridmap%area_src(ni)*re**2
     gdata_i = gdata_i + data_i(ni)*gridmap%area_src(ni)*gridmap%frac_src(ni)*re**2
  end do

  ! Sums on output grid

  gdata_o = 0.
  garea_o = 0.
  do no = 1,ns_o
     garea_o = garea_o + gridmap%area_dst(no)*re**2
     gdata_o = gdata_o + data_o(no)*gridmap%area_dst(no)*gridmap%frac_dst(no)*re**2
  end do

  ! Correct units

  if (percent) then
     gdata_i = gdata_i / 100._r8
     gdata_o = gdata_o / 100._r8
  end if

  ! Diagnostic output

  write (ndiag,*)
  write (ndiag,'(1x,70a1)') ('=',k=1,70)
  write (ndiag,*) trim(name), ' Output'
  write (ndiag,'(1x,70a1)') ('=',k=1,70)

  write (ndiag,*)
  write (ndiag,'(1x,70a1)') ('.',k=1,70)
  write (ndiag,2001)
2001 format (1x,'surface type   input grid area  output grid area'/ &
             1x,'                 10**6 km**2      10**6 km**2   ')
  write (ndiag,'(1x,70a1)') ('.',k=1,70)
  write (ndiag,*)
  write (ndiag,2002) name,          gdata_i*1.e-06, gdata_o*1.e-06
  write (ndiag,2002) 'all surface', garea_i*1.e-06, garea_o*1.e-06
2002 format (1x,a12,           f14.3,f17.3)

end subroutine output_diagnostics_area
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: output_diagnostics_continuous
!
! !INTERFACE:
subroutine output_diagnostics_continuous(data_i, data_o, gridmap, name, units, ndiag)
!
! !DESCRIPTION:
! Output diagnostics for a continuous field (but not area, for which there is a different routine)
!
! !USES:
  use mkgridmapMod, only : gridmap_type
  use mkvarpar, only : re
!
! !ARGUMENTS:
  implicit none
  real(r8)          , intent(in) :: data_i(:)    ! data on input grid
  real(r8)          , intent(in) :: data_o(:)    ! data on output grid
  type(gridmap_type), intent(in) :: gridmap      ! mapping info
  character(len=*)  , intent(in) :: name         ! name of field
  character(len=*)  , intent(in) :: units        ! units of field
  integer           , intent(in) :: ndiag        ! unit number for diagnostic output
!
! !REVISION HISTORY:
! Author: Bill Sacks
!
!
! !LOCAL VARIABLES:
!EOP
  real(r8) :: gdata_i         ! global sum of input data
  real(r8) :: gdata_o         ! global sum of output data
  real(r8) :: gwt_i           ! global sum of input weights (area * frac)
  real(r8) :: gwt_o           ! global sum of output weights (area * frac)
  integer  :: ns_i, ns_o      ! sizes of input & output grids
  integer  :: ni,no,k         ! indices

  character(len=*), parameter :: subname = "output_diagnostics_continuous"
!------------------------------------------------------------------------------

  ! Error check for array size consistencies

  ns_i = gridmap%na
  ns_o = gridmap%nb
  if (size(data_i) /= ns_i .or. &
      size(data_o) /= ns_o) then
     write(6,*) subname//' ERROR: array size inconsistencies for ', trim(name)
     write(6,*) 'size(data_i) = ', size(data_i)
     write(6,*) 'ns_i         = ', ns_i
     write(6,*) 'size(data_o) = ', size(data_o)
     write(6,*) 'ns_o         = ', ns_o
     stop
  end if

  ! Sums on input grid

  gdata_i = 0.
  gwt_i = 0.
  do ni = 1,ns_i
     gdata_i = gdata_i + data_i(ni)*gridmap%area_src(ni)*gridmap%frac_src(ni)
     gwt_i = gwt_i + gridmap%area_src(ni)*gridmap%frac_src(ni)
  end do

  ! Sums on output grid

  gdata_o = 0.
  gwt_o = 0.
  do no = 1,ns_o
     gdata_o = gdata_o + data_o(no)*gridmap%area_dst(no)*gridmap%frac_dst(no)
     gwt_o = gwt_o + gridmap%area_dst(no)*gridmap%frac_dst(no)
  end do

  ! Correct units

  gdata_i = gdata_i / gwt_i
  gdata_o = gdata_o / gwt_o

  ! Diagnostic output

  write (ndiag,*)
  write (ndiag,'(1x,70a1)') ('=',k=1,70)
  write (ndiag,*) trim(name), ' Output'
  write (ndiag,'(1x,70a1)') ('=',k=1,70)

  write (ndiag,*)
  write (ndiag,'(1x,70a1)') ('.',k=1,70)
  write (ndiag,2001)
  write (ndiag,2002) units, units
2001 format (1x,'   parameter              input grid          output grid')
2002 format (1x,'                 ', a24, a24)
  write (ndiag,'(1x,70a1)') ('.',k=1,70)
  write (ndiag,*)
  write (ndiag,2003) name,          gdata_i, gdata_o
2003 format (1x,a12,           f22.3,f17.3)

end subroutine output_diagnostics_continuous
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: output_diagnostics_continuous_outonly
!
! !INTERFACE:
subroutine output_diagnostics_continuous_outonly(data_o, gridmap, name, units, ndiag)
!
! !DESCRIPTION:
! Output diagnostics for a continuous field, just on the output grid
! This is used when the average of the field on the input grid is not of interest (e.g.,
! when the output quantity is the standard deviation of the input field)
!
! !USES:
  use mkgridmapMod, only : gridmap_type
  use mkvarpar, only : re
!
! !ARGUMENTS:
  implicit none
  real(r8)          , intent(in) :: data_o(:)    ! data on output grid
  type(gridmap_type), intent(in) :: gridmap      ! mapping info
  character(len=*)  , intent(in) :: name         ! name of field
  character(len=*)  , intent(in) :: units        ! units of field
  integer           , intent(in) :: ndiag        ! unit number for diagnostic output
!
! !REVISION HISTORY:
! Author: Bill Sacks
!
!
! !LOCAL VARIABLES:
!EOP
  real(r8) :: gdata_o         ! global sum of output data
  real(r8) :: gwt_o           ! global sum of output weights (area * frac)
  integer  :: ns_o            ! size of output grid
  integer  :: no,k            ! indices

  character(len=*), parameter :: subname = "output_diagnostics_continuous_outonly"
!------------------------------------------------------------------------------

  ! Error check for array size consistencies

  ns_o = gridmap%nb
  if (size(data_o) /= ns_o) then
     write(6,*) subname//' ERROR: array size inconsistencies for ', trim(name)
     write(6,*) 'size(data_o) = ', size(data_o)
     write(6,*) 'ns_o         = ', ns_o
     stop
  end if

  ! Sums on output grid

  gdata_o = 0.
  gwt_o = 0.
  do no = 1,ns_o
     gdata_o = gdata_o + data_o(no)*gridmap%area_dst(no)*gridmap%frac_dst(no)
     gwt_o = gwt_o + gridmap%area_dst(no)*gridmap%frac_dst(no)
  end do

  ! Correct units

  gdata_o = gdata_o / gwt_o

  ! Diagnostic output

  write (ndiag,*)
  write (ndiag,'(1x,70a1)') ('=',k=1,70)
  write (ndiag,*) trim(name), ' Output'
  write (ndiag,'(1x,70a1)') ('=',k=1,70)

  write (ndiag,*)
  write (ndiag,'(1x,70a1)') ('.',k=1,70)
  write (ndiag,2001)
  write (ndiag,2002) units
2001 format (1x,'   parameter              output grid')
2002 format (1x,'                 ', a24)
  write (ndiag,'(1x,70a1)') ('.',k=1,70)
  write (ndiag,*)
  write (ndiag,2003) name,          gdata_o
2003 format (1x,a12,           f22.3)

end subroutine output_diagnostics_continuous_outonly
!------------------------------------------------------------------------------


end module mkdiagnosticsMod
