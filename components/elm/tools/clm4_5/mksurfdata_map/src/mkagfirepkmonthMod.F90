module mkagfirepkmonthMod
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mkagfirepkmonthMod
!
! !DESCRIPTION:
! Make agricultural fire peak month data
!
! !REVISION HISTORY:
! Author: Sam Levis and Bill Sacks
!
!-----------------------------------------------------------------------
!!USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_flush
  use mkdomainMod , only : domain_checksame
  implicit none

  SAVE
  private           ! By default make data private
!
! !PUBLIC MEMBER FUNCTIONS:
!
  public mkagfirepkmon        ! Set agricultural fire peak month
!
! !PRIVATE MEMBER FUNCTIONS:
  private define_months       ! define month strings
!
! !PRIVATE DATA MEMBERS:
!
  integer , parameter :: min_valid_value = 1
  integer , parameter :: max_valid_value = 12
  integer , parameter :: unsetmon        = 13   ! flag to indicate agricultural fire peak month NOT set
!
! !PRIVATE DATA MEMBERS:
!
!EOP
!===============================================================
contains
!===============================================================

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkagfirepkmon
!
! !INTERFACE:
subroutine mkagfirepkmon(ldomain, mapfname, datfname, ndiag, &
                         agfirepkmon_o)
!
! !DESCRIPTION:
! Make agricultural fire peak month data from higher resolution data
!
! !USES:
  use mkdomainMod, only : domain_type, domain_clean, domain_read
  use mkgridmapMod
  use mkindexmapMod, only : get_dominant_indices
  use mkvarpar, only : re
  use mkncdio
  use mkchecksMod, only : min_bad, max_bad
!
! !ARGUMENTS:
  implicit none
  type(domain_type) , intent(in) :: ldomain
  character(len=*)  , intent(in) :: mapfname           ! input mapping file name
  character(len=*)  , intent(in) :: datfname           ! input data file name
  integer           , intent(in) :: ndiag              ! unit number for diag out
  integer           , intent(out):: agfirepkmon_o(:)   ! agricultural fire peak month
!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! Author: Sam Levis and Bill Sacks
!
!
! !LOCAL VARIABLES:
!EOP
  type(gridmap_type)    :: tgridmap
  type(domain_type)     :: tdomain          ! local domain
  real(r8), allocatable :: gast_i(:)        ! global area, by surface type
  real(r8), allocatable :: gast_o(:)        ! global area, by surface type
  integer , allocatable :: agfirepkmon_i(:) ! input grid: agricultural fire peak month
  integer  :: nagfirepkmon                  ! number of peak months 
  character(len=35), allocatable :: month(:)! name of each month
  integer  :: k,ni,no,ns_i,ns_o             ! indices
  integer  :: ncid,varid                    ! input netCDF id's
  integer  :: ier                           ! error status

  integer, parameter :: miss   = unsetmon   ! missing data indicator
  integer, parameter :: min_valid = 1       ! minimum valid value
  integer, parameter :: max_valid = 13      ! maximum valid value
  character(len=32) :: subname = 'mkagfirepkmon'
!-----------------------------------------------------------------------

  write (6,*) 'Attempting to make agricultural fire peak month data .....'
  call shr_sys_flush(6)

  ! -----------------------------------------------------------------
  ! Read domain and mapping information, check for consistency
  ! -----------------------------------------------------------------

  call domain_read( tdomain,datfname )

  call gridmap_mapread( tgridmap, mapfname )
  call gridmap_check( tgridmap, subname )

  call domain_checksame( tdomain, ldomain, tgridmap )

  ns_i = tdomain%ns
  ns_o = ldomain%ns

  ! -----------------------------------------------------------------
  ! Open input file, allocate memory for input data
  ! -----------------------------------------------------------------

  write (6,*) 'Open agricultural fire peak month file: ', trim(datfname)
  call check_ret(nf_open(datfname, 0, ncid), subname)

  allocate(agfirepkmon_i(ns_i), stat=ier)
  if (ier/=0) call abort()

  ! -----------------------------------------------------------------
  ! Regrid ag fire peak month
  ! -----------------------------------------------------------------

  call check_ret(nf_inq_varid (ncid, 'abm', varid), subname)
  call check_ret(nf_get_var_int (ncid, varid, agfirepkmon_i), subname)
  ! Note that any input point that is outside the range [min_valid_value,max_valid_value]
  ! will be ignored; this ignores input points with value of unsetmon
  call get_dominant_indices(tgridmap, agfirepkmon_i, agfirepkmon_o, &
       min_valid_value, max_valid_value, miss)

  ! Check validity of output data
  if (min_bad(agfirepkmon_o, min_valid, 'agfirepkmon') .or. &
      max_bad(agfirepkmon_o, max_valid, 'agfirepkmon')) then
     stop
  end if
  

  ! -----------------------------------------------------------------
  ! Output diagnostics comparing global area of each peak month on input and output grids
  !
  ! WJS (3-4-13): I am trying to generally put these diagnostics in mkdiagnosticsMod, but
  ! so far there isn't a general diagnostics routine for categorical data
  ! -----------------------------------------------------------------

  nagfirepkmon = maxval(agfirepkmon_i)
  allocate(gast_i(1:nagfirepkmon),gast_o(1:nagfirepkmon),month(1:nagfirepkmon))
  call define_months(nagfirepkmon, month)

  gast_i(:) = 0.0_r8
  do ni = 1,ns_i
     k = agfirepkmon_i(ni)
     gast_i(k) = gast_i(k) + tgridmap%area_src(ni)*tgridmap%frac_src(ni)*re**2
  end do

  gast_o(:) = 0.0_r8
  do no = 1,ns_o
     k = agfirepkmon_o(no)
     gast_o(k) = gast_o(k) + tgridmap%area_dst(no)*tgridmap%frac_dst(no)*re**2
  end do

  ! area comparison

  write (ndiag,*)
  write (ndiag,'(1x,70a1)') ('=',k=1,70)
  write (ndiag,*) 'Agricultural fire peak month Output'
  write (ndiag,'(1x,70a1)') ('=',k=1,70)

  write (ndiag,*)
  write (ndiag,'(1x,70a1)') ('.',k=1,70)
  write (ndiag,1001)
1001 format (1x,'peak month',20x,' input grid area output grid area',/ &
       1x,33x,'     10**6 km**2','      10**6 km**2')
  write (ndiag,'(1x,70a1)') ('.',k=1,70)
  write (ndiag,*)

  do k = 1, nagfirepkmon
     write (ndiag,1002) month(k),gast_i(k)*1.e-6,gast_o(k)*1.e-6
1002 format (1x,a35,f16.3,f17.3)
  end do

  ! -----------------------------------------------------------------
  ! Close files and deallocate dynamic memory
  ! -----------------------------------------------------------------

  call check_ret(nf_close(ncid), subname)
  call domain_clean(tdomain)
  call gridmap_clean(tgridmap)
  deallocate (agfirepkmon_i,gast_i,gast_o,month)

  write (6,*) 'Successfully made Agricultural fire peak month'
  write (6,*)
  call shr_sys_flush(6)

end subroutine mkagfirepkmon

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: define_months
!
! !INTERFACE:
subroutine define_months(nagfirepkmon, month)
!
! !DESCRIPTION:
! Define month strings
!
! !USES:
!
! !ARGUMENTS:
  implicit none
  integer         , intent(in) :: nagfirepkmon    ! max input value (including the 'unset' special value)
  character(len=*), intent(out):: month(:)        ! name of each month value
!
! !CALLED FROM:
! subroutine mkagfirepkmon
!
! !REVISION HISTORY:
! Author: Bill Sacks
!
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------

  if (nagfirepkmon == unsetmon) then
     if (size(month) < 13) then
        write(6,*) 'month array too small: ', size(month), ' < 13'
        call abort()
     end if
     month(1)  = 'January                             '
     month(2)  = 'February                            '
     month(3)  = 'March                               '
     month(4)  = 'April                               '
     month(5)  = 'May                                 '
     month(6)  = 'June                                '
     month(7)  = 'July                                '
     month(8)  = 'August                              '
     month(9)  = 'September                           '
     month(10) = 'October                             '
     month(11) = 'November                            '
     month(12) = 'December                            '
     month(13) = 'no agricultural fire peak month data'
  else
     write(6,*)'nagfirepkmon value of ',nagfirepkmon,' not supported'
     call abort()
  end if

end subroutine define_months
!-----------------------------------------------------------------------


end module mkagfirepkmonthMod
