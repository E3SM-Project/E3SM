module dynInitColumnsMod

  !---------------------------------------------------------------------------
  !
  ! !DESCRIPTION:
  ! Handle initialization of columns that just switched from inactive to active
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use decompMod         , only : bounds_type
  use abortutils        , only : endrun
  use clm_varctl        , only : iulog
  use elm_varcon        , only : ispval, namec
  use TemperatureType   , only : temperature_type
  use SoilHydrologyType , only : soilhydrology_type
  use WaterstateType    , only : waterstate_type
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
  subroutine initialize_new_columns(bounds, cactive_prior, &
       temperature_vars, waterstate_vars, soilhydrology_vars)
    !
    ! !DESCRIPTION:
    ! Do initialization for all columns that are newly-active in this time step
    !
    ! !USES:
    use GetGlobalValuesMod , only : GetGlobalWrite
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds                        ! bounds
    logical                  , intent(in)    :: cactive_prior( bounds%begc: ) ! column-level active flags from prior time step
    type(temperature_type)   , intent(inout) :: temperature_vars
    type(waterstate_type)    , intent(inout) :: waterstate_vars
    type(soilhydrology_type) , intent(inout) :: soilhydrology_vars
    !
    ! !LOCAL VARIABLES:
    integer :: c          ! column index
    integer :: c_template ! index of template column

    character(len=*), parameter :: subname = 'initialize_new_columns'
    !-----------------------------------------------------------------------
    
    SHR_ASSERT_ALL((ubound(cactive_prior) == (/bounds%endc/)), errMsg(__FILE__, __LINE__))

    do c = bounds%begc, bounds%endc
       ! If this column is newly-active, then we need to initialize it using the routines in this module
       if (col_pp%active(c) .and. .not. cactive_prior(c)) then
          c_template = initial_template_col_dispatcher(bounds, c, cactive_prior(bounds%begc:bounds%endc))
          if (c_template /= ispval) then
             call copy_state(c, c_template, temperature_vars, waterstate_vars, soilhydrology_vars)
          else
             write(iulog,*) subname// ' WARNING: No template column found to initialize newly-active column'
             write(iulog,*) '-- keeping the state that was already in memory, possibly from arbitrary initialization'
             call GetGlobalWrite(decomp_index=c, clmlevel=namec)
          end if
       end if
    end do

  end subroutine initialize_new_columns


  !-----------------------------------------------------------------------
  function initial_template_col_dispatcher(bounds, c_new, cactive_prior) result(c_template)
    !
    ! !DESCRIPTION:
    ! Find column to use as a template for the given column that has newly become active;
    ! this is a dispatcher that calls the appropriate routine based on the landunit type of c_new.
    !
    ! Returns ispval if there is no column to use for initialization
    !
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

    character(len=*), parameter :: subname = 'initial_template_col_dispatcher'
    !-----------------------------------------------------------------------
    
    SHR_ASSERT_ALL((ubound(cactive_prior) == (/bounds%endc/)), errMsg(__FILE__, __LINE__))

    l = col_pp%landunit(c_new)
    ltype = lun_pp%itype(l)
    select case(ltype)
    case(istsoil)
       c_template = initial_template_col_soil(c_new)
    case(istcrop)
       c_template = initial_template_col_crop(bounds, c_new, cactive_prior(bounds%begc:bounds%endc))
    case(istice)
       write(iulog,*) subname// ' ERROR: Ability to initialize a newly-active glacier column not yet implemented'
       write(iulog,*) 'Expectation is that only ice_mec columns can grow'
       call endrun(decomp_index=c_new, clmlevel=namec, msg=errMsg(__FILE__, __LINE__))
    case(istice_mec)
       write(iulog,*) subname// ' ERROR: Ability to initialize a newly-active glacier mec column not yet implemented'
       write(iulog,*) 'Expectation is that glacier mec columns should be active from the start of the run wherever they can grow'
       call endrun(decomp_index=c_new, clmlevel=namec, msg=errMsg(__FILE__, __LINE__))
    case(istdlak)
       write(iulog,*) subname// ' ERROR: Ability to initialize a newly-active lake column not yet implemented'
       call endrun(decomp_index=c_new, clmlevel=namec, msg=errMsg(__FILE__, __LINE__))
    case(istwet)
       write(iulog,*) subname// ' ERROR: Ability to initialize a newly-active wetland column not yet implemented'
       call endrun(decomp_index=c_new, clmlevel=namec, msg=errMsg(__FILE__, __LINE__))
    case(isturb_MIN:isturb_MAX)
       write(iulog,*) subname// ' ERROR: Ability to initialize a newly-active urban column not yet implemented'
       call endrun(decomp_index=c_new, clmlevel=namec, msg=errMsg(__FILE__, __LINE__))
    case default
       write(iulog,*) subname// ' ERROR: Unknown landunit type: ', ltype
       call endrun(decomp_index=c_new, clmlevel=namec, msg=errMsg(__FILE__, __LINE__))
    end select

  end function initial_template_col_dispatcher


  !-----------------------------------------------------------------------
  function initial_template_col_soil(c_new) result(c_template)
    !
    ! !DESCRIPTION:
    ! Find column to use as a template for a vegetated column that has newly become active.
    !
    ! For now, we assume that the only vegetated columns that can newly become active are
    ! ones with 0 weight on the grid cell (i.e., virtual columns). For these, we simply
    ! keep the state at the current value (likely arbitrary initial conditions), and so
    ! return ispval from this function. Within this function, we check this assumption.
    !
    ! !USES:
    use elm_varcon, only : ispval
    !
    ! !ARGUMENTS:
    integer              :: c_template ! function result
    integer , intent(in) :: c_new        ! column index that needs initialization
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'initial_template_col_soil'
    !-----------------------------------------------------------------------

    if (col_pp%wtgcell(c_new) > 0._r8) then
       write(iulog,*) subname// ' ERROR: Expectation is that the only vegetated columns that&
            & can newly become active are ones with 0 weight on the grid cell'
       call endrun(decomp_index=c_new, clmlevel=namec, msg=errMsg(__FILE__, __LINE__))
    end if

    c_template = ispval
    
  end function initial_template_col_soil

  !-----------------------------------------------------------------------
  function initial_template_col_crop(bounds, c_new, cactive_prior) result(c_template)
    !
    ! !DESCRIPTION:
    ! Find column to use as a template for a crop column that has newly become active
    !
    ! Returns ispval if there is no column to use for initialization
    !
    ! !USES:
    use elm_varcon, only : ispval
    use landunit_varcon, only : istsoil, istcrop
    !
    ! !ARGUMENTS:
    integer :: c_template  ! function result
    type(bounds_type) , intent(in) :: bounds                        ! bounds
    integer           , intent(in) :: c_new                         ! column index that needs initialization
    logical           , intent(in) :: cactive_prior( bounds%begc: ) ! column-level active flags from prior time step
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'initial_template_col_crop'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(cactive_prior) == (/bounds%endc/)), errMsg(__FILE__, __LINE__))
    
    ! First try to find an active column on the vegetated landunit; if there is none, then
    ! find the first active column on the crop landunit; if there is none, then
    ! template_col will be ispval
    c_template = initial_template_col(bounds, c_new, istsoil, cactive_prior(bounds%begc:bounds%endc))
    if (c_template == ispval) then
       c_template = initial_template_col(bounds, c_new, istcrop, cactive_prior(bounds%begc:bounds%endc))
    end if

  end function initial_template_col_crop


  !-----------------------------------------------------------------------
  function initial_template_col(bounds, c_new, landunit_type, cactive_prior) result(c_template)
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
    !
    ! !USES:
    use elm_varcon, only : ispval
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
    
    character(len=*), parameter :: subname = 'initial_template_col'
    !-----------------------------------------------------------------------
    
    SHR_ASSERT_ALL((ubound(cactive_prior) == (/bounds%endc/)), errMsg(__FILE__, __LINE__))

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
  subroutine copy_state(c_new, c_template, &
       temperature_vars, waterstate_vars, soilhydrology_vars)
    !
    ! !DESCRIPTION:
    ! Copy a subset of state variables from a template column (c_template) to a newly-
    ! active column (c_new)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer                  , intent(in)    :: c_new      ! index of newly-active column
    integer                  , intent(in)    :: c_template ! index of column to use as a template
    type(temperature_type)   , intent(inout) :: temperature_vars
    type(waterstate_type)    , intent(inout) :: waterstate_vars
    type(soilhydrology_type) , intent(inout) :: soilhydrology_vars
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'copy_state'
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
