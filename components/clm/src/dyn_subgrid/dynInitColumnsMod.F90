module dynInitColumnsMod

  !---------------------------------------------------------------------------
  !
  ! !DESCRIPTION:
  ! Handle initialization of columns that just switched from inactive to active
  !
  ! !USES:
!#py #include "shr_assert.h"
  use shr_kind_mod      , only : r8 => shr_kind_r8
  !#py !#py use shr_log_mod       , only : errMsg => shr_log_errMsg
  use decompMod         , only : bounds_type
  !#py use abortutils        , only : endrun
  use clm_varctl        , only : iulog
  use clm_varcon        , only : ispval, namec
  use SoilHydrologyType , only : soilhydrology_type
  use TopounitType      , only : top_pp
  use LandunitType      , only : lun_pp
  use ColumnType        , only : col_pp
  use ColumnDataType    , only : col_es, col_ws
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  private
  save
  ! The following is the public interface to the routines in this module:
  public :: initialize_new_columns  ! Do initialization for all columns that are newly-active in this time step

  ! The following are public only for unit testing purposes, and should not be called
  ! directly by application code:
  public :: initial_template_col_crop ! Find column to use as a template for a crop column that has newly become active
  public :: initial_template_col      ! Find column to serve as a template for the initialization of another column in the same grid cell
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: initial_template_col_dispatcher ! Find column to use as a template; dispatcher to the appropriate routine based on landunit type
  private :: initial_template_col_soil       ! Find column to use as a template for a vegetated column that has newly become active
  private :: copy_state                      ! Copy a subset of state variables from template column to newly-active column

  !---------------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine initialize_new_columns(bounds, cactive_prior, soilhydrology_vars)
    !
    ! !DESCRIPTION:
    ! Do initialization for all columns that are newly-active in this time step
    !
    ! !USES:
    !
    ! !ARGUMENTS:
      !$acc routine seq
    type(bounds_type)        , intent(in)    :: bounds                        ! bounds
    logical                  , intent(in)    :: cactive_prior( bounds%begc: ) ! column-level active flags from prior time step
    type(soilhydrology_type) , intent(inout) :: soilhydrology_vars
    !
    ! !LOCAL VARIABLES:
    integer :: c          ! column index
    integer :: c_template ! index of template column

    !-----------------------------------------------------------------------


    do c = bounds%begc, bounds%endc
       ! If this column is newly-active, then we need to initialize it using the routines in this module
       if (col_pp%active(c) .and. .not. cactive_prior(c)) then
          c_template = initial_template_col_dispatcher(bounds, c, cactive_prior(bounds%begc:bounds%endc))
          if (c_template /= ispval) then
             call copy_state(c, c_template, soilhydrology_vars)
          else
             print *, ' WARNING: No template column found to initialize newly-active column'
          end if
       end if
    end do

  end subroutine initialize_new_columns


  !-----------------------------------------------------------------------
  integer function initial_template_col_dispatcher(bounds, c_new, cactive_prior) result(c_template)
    !
    ! !DESCRIPTION:
    ! Find column to use as a template for the given column that has newly become active;
    ! this is a dispatcher that calls the appropriate routine based on the landunit type of c_new.
    !
    ! Returns ispval if there is no column to use for initialization
    !$acc routine seq
    ! !USES:
    use landunit_varcon, only : istsoil, istcrop, istice, istice_mec, istdlak, istwet, isturb_MIN, isturb_MAX
    !
    ! !ARGUMENTS:
    integer :: c_template  ! function result
    type(bounds_type) , intent(in) :: bounds                        ! bounds
    integer           , intent(in) :: c_new                         ! column index that needs initialization
    logical           , intent(in) :: cactive_prior( bounds%begc: ) ! column-level active flags from prior time step
    !
    ! !LOCAL VARIABLES:
    integer :: l     ! landunit index
    integer :: ltype ! landunit type

    !-----------------------------------------------------------------------


    l = col_pp%landunit(c_new)
    ltype = lun_pp%itype(l)
    select case(ltype)
    case(istsoil)
       c_template = initial_template_col_soil(c_new)
    case(istcrop)
       c_template = initial_template_col_crop(bounds, c_new,cactive_prior(bounds%begc:bounds%endc) )
    case(istice)
       print *, ' ERROR: Ability to initialize a newly-active glacier column not yet implemented'
    case(istice_mec)
       print *, ' ERROR: Ability to initialize a newly-active glacier mec column not yet implemented'
    case(istdlak)
       print *, ' ERROR: Ability to initialize a newly-active lake column not yet implemented'
    case(istwet)
       print *, ' ERROR: Ability to initialize a newly-active wetland column not yet implemented'
    case(isturb_MIN:isturb_MAX)
       print *, ' ERROR: Ability to initialize a newly-active urban column not yet implemented'
    case default
       print *, ' ERROR: Unknown landunit type: ', ltype
    end select

  end function initial_template_col_dispatcher


  !-----------------------------------------------------------------------
  integer function initial_template_col_soil(c_new) result(c_template)
    !
    ! !DESCRIPTION:
    ! Find column to use as a template for a vegetated column that has newly become active.
    !
    ! For now, we assume that the only vegetated columns that can newly become active are
    ! ones with 0 weight on the grid cell (i.e., virtual columns). For these, we simply
    ! keep the state at the current value (likely arbitrary initial conditions), and so
    ! return ispval from this function. Within this function, we check this assumption.
    !$acc routine seq
    ! !USES:
    use clm_varcon, only : ispval
    !
    ! !ARGUMENTS:
    integer              :: c_template ! function result
    integer , intent(in) :: c_new        ! column index that needs initialization
    !
    ! !LOCAL VARIABLES:
    !-----------------------------------------------------------------------

    if (col_pp%wtgcell(c_new) > 0._r8) then

       print *, ' ERROR: Expectation is that the only vegetated columns that',c_new
       print *, ' can newly become active are ones with 0 weight on the grid cell'
       stop 
    end if

    c_template = ispval

  end function initial_template_col_soil

  !-----------------------------------------------------------------------
  integer function initial_template_col_crop(bounds, c_new, cactive_prior) result(c_template)
    !
    ! !DESCRIPTION:
    ! Find column to use as a template for a crop column that has newly become active
    !
    ! Returns ispval if there is no column to use for initialization
    !$acc routine seq
    ! !USES:
    use clm_varcon     , only : ispval
    use landunit_varcon, only : istsoil, istcrop
    !
    ! !ARGUMENTS:
    integer :: c_template  ! function result
    type(bounds_type) , intent(in) :: bounds                        ! bounds
    integer           , intent(in) :: c_new                         ! column index that needs initialization
    logical           , intent(in) :: cactive_prior( bounds%begc: ) ! column-level active flags from prior time step
    !
    ! !LOCAL VARIABLES:

    !-----------------------------------------------------------------------

    ! First try to find an active column on the vegetated landunit; if there is none, then
    ! find the first active column on the crop landunit; if there is none, then
    ! template_col will be ispval
    c_template = initial_template_col(bounds, c_new, istsoil, cactive_prior(bounds%begc:bounds%endc) )
    if (c_template == ispval) then
       c_template = initial_template_col(bounds, c_new, istcrop, cactive_prior(bounds%begc:bounds%endc) )
    end if

  end function initial_template_col_crop


  !-----------------------------------------------------------------------
  integer function initial_template_col(bounds, c_new, landunit_type, cactive_prior) result(c_template)
    !
    ! !DESCRIPTION:
    ! Finds a column to serve as a template for the initialization of another column in
    ! the same grid cell.
    !
    ! Looks for a landunit of the type given by landunit_type (e.g., istsoil,
    ! istcrop). Looks for the first active column on this landunit type, in the same grid
    ! cell; order of columns within a landunit is arbitrary (given by their order in
    ! memory). Returns the column index of the first such column found. If there are no
    ! active columns in this landunit in this grid cell, returns ispval.
    !
    ! Note that, in checking 'active', we use the active flags from the prior time step,
    ! so that we don't identify a point that just became active for the first time in this
    ! time step.
    !$acc routine seq
    ! !USES:
    use clm_varcon, only : ispval
    !
    ! !ARGUMENTS:
    integer :: c_template  ! function return value

    type(bounds_type) , intent(in) :: bounds                        ! bounds
    integer           , intent(in) :: c_new                         ! column index that needs initialization
    integer           , intent(in) :: landunit_type                 ! landunit type from which we want to find a template column (e.g., istsoil)
    logical           , intent(in) :: cactive_prior( bounds%begc: ) ! column-level active flags from prior time step
    !
    ! !LOCAL VARIABLES:
    logical :: found  ! whether a suitable template column has been found
    integer :: t,l,c  ! indices of topounit, landunit, column

    !-----------------------------------------------------------------------


    found = .false.
    t = col_pp%topounit(c_new)
    l = top_pp%landunit_indices(landunit_type, t)

    ! If this landunit exists on this grid cell...
    if (l /= ispval) then

       ! Loop through columns on this landunit; stop if as soon as we find an active
       ! column: that will serve as the template
       c = lun_pp%coli(l)
       do while (.not. found .and. c <= lun_pp%colf(l))
          if (cactive_prior(c)) then
             found = .true.
          else
             c = c + 1
          end if
       end do
    end if

    if (found) then
       c_template = c
    else
       c_template = ispval
    end if

  end function initial_template_col

  !-----------------------------------------------------------------------
  subroutine copy_state(c_new, c_template, soilhydrology_vars)
    !
    ! !DESCRIPTION:
    ! Copy a subset of state variables from a template column (c_template) to a newly-
    ! active column (c_new)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
      !$acc routine seq
    integer                  , intent(in)    :: c_new      ! index of newly-active column
    integer                  , intent(in)    :: c_template ! index of column to use as a template
    type(soilhydrology_type) , intent(inout) :: soilhydrology_vars
    !
    ! !LOCAL VARIABLES:

    !-----------------------------------------------------------------------

    ! For now, just copy t_soisno
    ! TODO: Figure out what else should be copied
    col_es%t_soisno(c_new,:) = col_es%t_soisno(c_template,:)

    ! TODO(wjs, 2016-08-31) If we had more general uses of this initial template col
    ! infrastructure (copying state between very different landunits), then we might need
    ! to handle bedrock layers - e.g., zeroing out any water that would be added to a
    ! bedrock layer(?). But for now we just use this initial template col infrastructure
    ! for nat veg -> crop, for which the bedrock will be the same, so we're not dealing
    ! with that complexity for now.
    col_ws%h2osoi_liq(c_new,1:) = col_ws%h2osoi_liq(c_template,1:)
    col_ws%h2osoi_ice(c_new,1:) = col_ws%h2osoi_ice(c_template,1:)
    col_ws%h2osoi_vol(c_new,1:) = col_ws%h2osoi_vol(c_template,1:)

    soilhydrology_vars%wa_col(c_new) = soilhydrology_vars%wa_col(c_template)

   end subroutine copy_state



end module dynInitColumnsMod
