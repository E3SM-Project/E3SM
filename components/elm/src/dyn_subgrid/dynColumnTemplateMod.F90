module dynColumnTemplateMod

  !---------------------------------------------------------------------------
  !
  ! !DESCRIPTION:
  ! Routines for finding a template column to use for the state variables on some other
  ! column of interest.
  !
  ! For example, if a glacier column (with no carbon/nitrogen information) shrinks, and we
  ! want to assume that its carbon and nitrogen state implicitly matches the state of the
  ! vegetated landunit in that grid cell, then we can use a routine in this module to find
  ! that vegetated lanadunit.
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use decompMod       , only : bounds_type
  use GridcellType    , only : grc_pp
  use LandunitType    , only : lun_pp
  use ColumnType      , only : col_pp
  use elm_varcon      , only : ispval

  !
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  private

  ! ------------------------------------------------------------------------
  ! Functions that operate on a single column at a time
  ! ------------------------------------------------------------------------

  ! Find column to use as template by looking for an active column on a particular landunit
  public :: template_col_from_landunit

  ! ------------------------------------------------------------------------
  ! Subroutines that operate on the whole column-level array at once
  ! ------------------------------------------------------------------------

  ! Find column to use as template by looking for an active column on the natural veg landunit
  public :: template_col_from_natveg_array

  !
  ! !PUBLIC VARIABLES:

  ! if no template column was found, this value is returned
  integer, parameter, public :: TEMPLATE_NONE_FOUND = ispval

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  ! ------------------------------------------------------------------------
  ! Functions that operate on a single column at a time
  ! ------------------------------------------------------------------------

  !-----------------------------------------------------------------------
  function template_col_from_landunit(bounds, c_target, landunit_type, cactive) result(c_template)
    !
    ! !DESCRIPTION:
    ! Finds a column to serve as a template for the state variables on the target column.
    !
    ! Looks for a landunit of the type given by landunit_type (e.g., istsoil,
    ! istcrop). Looks for the first active column on this landunit type, in the same grid
    ! cell; order of columns within a landunit is arbitrary (given by their order in
    ! memory). Returns the column index of the first such column found. If there are no
    ! active columns in this landunit in this grid cell, returns TEMPLATE_NONE_FOUND.
    !
    ! Note that it is often most appropriate for cactive to be the active flags from the
    ! *prior* time step, so that we don't identify a point that just became active for the
    ! first time in this time step.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer :: c_template  ! function return value

    type(bounds_type) , intent(in) :: bounds                  ! bounds
    integer           , intent(in) :: c_target                ! column index for which we want a template
    integer           , intent(in) :: landunit_type           ! landunit type from which we want to find a template column (e.g., istsoil)
    logical           , intent(in) :: cactive( bounds%begc: ) ! column-level active flags (generally from prior time step)
    !
    ! !LOCAL VARIABLES:
    logical :: found  ! whether a suitable template column has been found
    integer :: g,l,c  ! indices of grid cell, landunit, column
    
    character(len=*), parameter :: subname = 'template_col_from_landunit'
    !-----------------------------------------------------------------------
    
    SHR_ASSERT_ALL((ubound(cactive) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))

    found = .false.
    g = col_pp%gridcell(c_target)
    l = grc_pp%landunit_indices(landunit_type, g)

    ! If this landunit exists on this grid cell...
    if (l /= ispval) then

       ! Loop through columns on this landunit; stop if as soon as we find an active
       ! column: that will serve as the template
       c = lun_pp%coli(l)
       do while (.not. found .and. c <= lun_pp%colf(l))
          if (cactive(c)) then
             found = .true.
          else
             c = c + 1
          end if
       end do
    end if

    if (found) then
       c_template = c
    else
       c_template = TEMPLATE_NONE_FOUND
    end if

  end function template_col_from_landunit

  ! ------------------------------------------------------------------------
  ! Subroutines that operate on the whole column-level array at once
  ! ------------------------------------------------------------------------

  !-----------------------------------------------------------------------
  subroutine template_col_from_natveg_array(bounds, cactive, c_templates)
    !
    ! !DESCRIPTION:
    ! For each column, finds a column to serve as a template for the state variables on
    ! the target column.
    !
    ! For each column: Looks for the first active column on the natural veg landunit in
    ! the same grid cell as the target column. If there are no active columns in the
    ! natural veg landunit in this grid cell, assigns TEMPLATE_NONE_FOUND.
    !
    ! Note: If there are multiple columns on the natural veg. landunit, then a given
    ! natural veg column may have a template col that differs from itself! The caller is
    ! responsible for determining when this template column should be used and when it
    ! should not be used.
    !
    ! See also the notes about cactive under 'template_col_from_landunit'.
    !
    ! !USES:
    use landunit_varcon, only : istsoil
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in)  :: bounds                      ! bounds
    logical           , intent(in)  :: cactive( bounds%begc: )     ! column-level active flags (generally from prior time step)
    integer           , intent(out) :: c_templates( bounds%begc: ) ! template column for each column
    !
    ! !LOCAL VARIABLES:
    integer :: c

    character(len=*), parameter :: subname = 'template_col_from_natveg_array'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(cactive) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(c_templates) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))

    do c = bounds%begc, bounds%endc
       c_templates(c) = template_col_from_landunit(bounds, c, istsoil, &
            cactive(bounds%begc:bounds%endc))
    end do

  end subroutine template_col_from_natveg_array

  

end module dynColumnTemplateMod
