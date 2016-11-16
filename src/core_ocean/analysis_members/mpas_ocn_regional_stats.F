! Copyright (c) 2015,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! ocn_regional_stats
!
!> \brief MPAS ocean analysis core member: regional_stats
!> \author Jon Woodring
!> \date   December 17, 2015
!> \details
!>  Flexible regional averaging, mins, and maxes of fields.
!-----------------------------------------------------------------------
module ocn_regional_stats
  use mpas_derived_types
  use mpas_pool_routines
  use mpas_dmpar
  use mpas_timekeeping
  use mpas_stream_manager

  use ocn_constants
  use ocn_diagnostics_routines

  implicit none
  private
  save

  ! Public parameters
  !--------------------------------------------------------------------

  ! Public member functions
  !--------------------------------------------------------------------
  public :: &
         ocn_init_regional_stats, &
         ocn_compute_regional_stats, &
         ocn_restart_regional_stats, &
         ocn_finalize_regional_stats

  ! Private module variables
  !--------------------------------------------------------------------

  type regional_variable_type
    ! state per instance variable, stored in framework
    character (len=StrKIND), pointer :: input_name
    integer, pointer :: has_vertical

    ! generated on every instance call
    character (len=StrKIND), dimension(:), allocatable :: output_names
  end type regional_variable_type

  type regional_type
    ! state per instance, stored in framework
    integer, pointer :: number_of_variables

    integer, pointer :: operation
    integer, pointer :: region_element
    integer, pointer :: function_oned
    integer, pointer :: function_twod
    integer, pointer :: group_index

    ! looked up on every instance call
    character (len=StrKIND), pointer :: weights_oned
    character (len=StrKIND), pointer :: weights_twod
    integer, dimension(:), pointer :: num_regions_per
    integer, dimension(:, :), pointer :: groups
    character (len=StrKIND), pointer :: vertical_mask, vertical_dim

    ! generated on every instance call
    character (len=StrKIND) :: masking_field

    ! allocated on every instance call
    type (regional_variable_type), dimension(:), allocatable :: variables
    character (len=StrKIND), dimension(:), allocatable :: count_zerod_names
    character (len=StrKIND), dimension(:), allocatable :: weight_zerod_names
    character (len=StrKIND), dimension(:), allocatable :: count_oned_names
    character (len=StrKIND), dimension(:), allocatable :: weight_oned_names
  end type regional_type

  ! enum of ops and types
  integer, parameter :: AVG_OP = 1
  integer, parameter :: MIN_OP = 2
  integer, parameter :: MAX_OP = 3
  integer, parameter :: SUM_OP = 4
  integer, parameter :: SOS_OP = 5

  integer, parameter :: ID_FUNC = 11
  integer, parameter :: MUL_FUNC = 12

  integer, parameter :: CELL_REGION = 101
  integer, parameter :: VERTEX_REGION = 102

  ! namelist operators and identifiers for unique names
  character (len=3), parameter :: AVG_TOKEN = 'avg'
  character (len=3), parameter :: MIN_TOKEN = 'min'
  character (len=3), parameter :: MAX_TOKEN = 'max'
  character (len=3), parameter :: SUM_TOKEN = 'sum'
  character (len=3), parameter :: SOS_TOKEN = 'sos'

  character (len=StrKIND), parameter :: COUNT_TOKEN = 'count'
  character (len=StrKIND), parameter :: WEIGHT_TOKEN = 'weight'
  character (len=StrKIND), parameter :: ZEROD_TOKEN = '0D'
  character (len=StrKIND), parameter :: ONED_TOKEN = '1D'

  ! namelist operators
  character (len=2), parameter :: ID_TOKEN  = 'id'
  character (len=3), parameter :: MUL_TOKEN = 'mul'

  ! namelist operators and identifiers for mask-struct data
  character (len=5), parameter :: CELL_TOKEN = 'cell'
  character (len=8), parameter :: VERTEX_TOKEN = 'vertex'

  ! camel case token for mask field names
  character (len=5), parameter :: CELL_CAMEL = 'Cell'
  character (len=8), parameter :: VERTEX_CAMEL = 'Vertex'

  ! canonical dimension names
  character (len=6), parameter :: CELL_DIM = 'nCells'
  character (len=9), parameter :: VERTEX_DIM = 'nVertices'

  ! canonical solve names
  character (len=11), parameter :: CELL_SOLVE = 'nCellsSolve'
  character (len=14), parameter :: VERTEX_SOLVE = 'nVerticesSolve'

  ! canonical mins and maxes
  real (kind=RKIND), parameter :: DEFAULT_MPAS_MIN_VALUE = -1.0e34_RKIND
  real (kind=RKIND), parameter :: DEFAULT_MPAS_MAX_VALUE =  1.0e34_RKIND

  ! none
  character (len=4), parameter :: NONE_TOKEN = 'none'

  ! memory names for duplication
  character (len=StrKIND), parameter :: REGIONAL_STATS_POOL = &
    'regionalStatsAM'
  character (len=StrKIND), parameter :: ONE_INTEGER_MEMORY = &
    'regionalStatsOneInteger'
  character (len=StrKIND), parameter :: ONE_STRING_MEMORY = &
    'regionalStatsOneString'
  character (len=StrKIND), parameter :: ONE_REAL_MEMORY = &
    'regionalStatsOneReal'

  ! mask-struct array names and suffixes
  character (len=StrKIND), parameter :: MASK_POOL_NAME = 'regions'
  character (len=StrKIND), parameter :: MASK_DATA_PREFIX = 'region'
  character (len=StrKIND), parameter :: MASK_DATA_SUFFIX = 'Masks'
  character (len=StrKIND), parameter :: GROUP_DATA_NAME = &
    'regionsInGroup'
  character (len=StrKIND), parameter :: REGIONS_PER_NAME = &
    'nRegionsInGroup'
  character (len=StrKIND), parameter :: REGION_NAMES_NAME = &
    'regionNames'
  character (len=StrKIND), parameter :: GROUP_NAMES_NAME = &
    'regionGroupNames'

  ! mask dimension names
  character (len=StrKIND), parameter :: NUM_REGIONS_SUFFIX = 'nRegions'
  character (len=StrKIND), parameter :: NUM_GROUPS_SUFFIX = 'nRegionGroups'
  character (len=StrKIND), parameter :: MAX_REGIONS_SUFFIX = &
    'maxRegionsInGroup'

  ! prefixes
  character (len=StrKIND), parameter :: CONFIG_PREFIX = &
    'config_AM_regionalStats'
  character (len=StrKIND), parameter :: FRAMEWORK_PREFIX = 'regionalStats'

  ! namelist-only suffixes
  character (len=StrKIND), parameter :: OUTPUT_STREAM_SUFFIX = '_output_stream'
  character (len=StrKIND), parameter :: WEIGHTS_ONED_SUFFIX = &
    '_1d_weighting_field'
  character (len=StrKIND), parameter :: WEIGHTS_TWOD_SUFFIX = &
    '_2d_weighting_field'
  character (len=StrKIND), parameter :: REGION_GROUP_SUFFIX = '_region_group'
  character (len=StrKIND), parameter :: INPUT_STREAM_SUFFIX = '_input_stream'
  character (len=StrKIND), parameter :: RESTART_STREAM_SUFFIX = &
    '_restart_stream'
  character (len=StrKIND), parameter :: VERTICAL_MASK_SUFFIX = &
    '_vertical_mask'
  character (len=StrKIND), parameter :: VERTICAL_DIM_SUFFIX = &
    '_vertical_dimension'

  ! namelist and instance suffixes
  character (len=StrKIND), parameter :: OPERATION_SUFFIX = '_operation'
  character (len=StrKIND), parameter :: REGION_TYPE_SUFFIX = '_region_type'
  character (len=StrKIND), parameter :: FUNCTION_ONED_SUFFIX = &
    '_1d_weighting_function'
  character (len=StrKIND), parameter :: FUNCTION_TWOD_SUFFIX = &
    '_2d_weighting_function'

  ! instance-only suffixes
  character (len=StrKIND), parameter :: INPUT_NAME_SUFFIX = '_input_name'
  character (len=StrKIND), parameter :: NUMBER_OF_VARIABLES_SUFFIX = &
    '_number_of_variables'
  character (len=StrKIND), parameter :: HAS_VERTICAL_SUFFIX = '_has_vertical'

  ! error message
  character (len=StrKIND), parameter :: CURRENT_CORE_NAME = 'MPAS-Ocean'

!***********************************************************************
contains



!***********************************************************************
! routine ocn_init_regional_stats
!
!> \brief Initialize MPAS-Ocean analysis member
!> \author  Jon Woodring
!> \date    December 17, 2015
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Ocean analysis member.
!-----------------------------------------------------------------------
subroutine ocn_init_regional_stats(domain, instance, err)!{{{
  ! input variables
  character (len=StrKIND), intent(in) :: instance

  ! input/output variables
  type (domain_type), intent(inout) :: domain

  ! output variables
  integer, intent(out) :: err !< Output: error flag

  ! local variables
  integer :: v

  ! start procedure
  err = 0

  ! create all of the state for this instance
  call start_state(domain, instance, err)

end subroutine ocn_init_regional_stats!}}}



!***********************************************************************
! routine ocn_compute_regional_stats
!
!> \brief Compute MPAS-Ocean analysis member
!> \author  Jon Woodring
!> \date    December 17, 2015
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!-----------------------------------------------------------------------
subroutine ocn_compute_regional_stats(domain, timeLevel, instance, err)!{{{
  ! input variables
  integer, intent(in) :: timeLevel
  character (len=StrKIND), intent(in) :: instance

  ! input/output variables
  type (domain_type), intent(inout) :: domain

  ! output variables
  integer, intent(out) :: err !< Output: error flag

  ! local variables
  integer :: v, last, m, b, i
  integer, pointer :: levels, solve
  type (regional_type) :: regions
  integer, pointer :: count_sum_zerod
  integer, dimension(:), pointer :: count_sum_oned
  real (kind=RKIND), pointer :: weight_sum_zerod
  real (kind=RKIND), dimension(:), pointer :: weight_sum_oned
  type (block_type), pointer :: block
  integer, dimension(:,:), pointer :: mask, vertical_mask
  type (mpas_pool_type), pointer :: amPool, maskPool
  real (kind=RKIND), dimension(:), pointer :: oned_weights
  real (kind=RKIND), dimension(:,:), pointer :: twod_weights
  logical :: active_vertical
  real (kind=RKIND) :: real_copy
  integer :: int_copy
  real (kind=RKIND), dimension(:), allocatable :: real_vert_copy
  integer, dimension(:), allocatable :: int_vert_copy

  ! start procedure
  err = 0

  ! get all of the state for this instance to be able to compute
  call get_state(domain, instance, regions)

  !
  ! calculate all of the counts and weights at the beginning
  call mpas_pool_get_subpool(domain % blocklist % structs, &
    REGIONAL_STATS_POOL, amPool)

  active_vertical = (trim(regions % vertical_dim) /= trim(NONE_TOKEN))
  if (active_vertical) then
    call mpas_pool_get_dimension(domain % blocklist % dimensions, &
      regions % vertical_dim, levels)
    allocate(real_vert_copy(levels))
    allocate(int_vert_copy(levels))
  end if

  last = regions % num_regions_per(regions % group_index)
  do b = 1, last
    m = regions % groups(b, regions % group_index)

    ! get the different target arrays
    ! count 0d
    call mpas_pool_get_array(amPool, &
      regions % count_zerod_names(b), count_sum_zerod, 1)
    count_sum_zerod = 0

    ! weight 0d
    if (regions % function_oned == MUL_FUNC) then
      call mpas_pool_get_array(amPool, &
        regions % weight_zerod_names(b), weight_sum_zerod, 1)
      weight_sum_zerod = 0
    end if

    ! count 1d
    if (active_vertical) then
      call mpas_pool_get_array(amPool, &
        regions % count_oned_names(b), count_sum_oned, 1)
      count_sum_oned = 0

      ! weight 1d
      if (regions % function_twod == MUL_FUNC) then
        call mpas_pool_get_array(amPool, &
          regions % weight_oned_names(b), weight_sum_oned, 1)
        weight_sum_oned = 0
      end if
    end if

    ! iterate over all the blocks and sum mask counts and weights
    block => domain % blocklist
    do while (associated(block))
      ! get the mask
      call mpas_pool_get_subpool(block % structs, MASK_POOL_NAME, maskPool)
      call mpas_pool_get_array(maskPool, regions % masking_field, mask, 1)

      ! get the dimension
      if (regions % region_element == CELL_REGION) then
        call mpas_pool_get_dimension(block % dimensions, CELL_SOLVE, solve)
      else
        call mpas_pool_get_dimension(block % dimensions, VERTEX_SOLVE, solve)
      end if

      ! create the counts and weights
      do i = 1, solve
        count_sum_zerod = count_sum_zerod + mask(m, i)
      end do

      if (regions % function_oned == MUL_FUNC) then
        call mpas_pool_get_array(block % allFields, &
          regions % weights_oned, oned_weights, 1)

        do i = 1, solve
          weight_sum_zerod = weight_sum_zerod + mask(m, i) * oned_weights(i)
        end do
      end if

      if (active_vertical) then
        call mpas_pool_get_array(block % allFields, &
          regions % vertical_mask, vertical_mask, 1)

        do i = 1, solve
          do v = 1, levels
            count_sum_oned(v) = count_sum_oned(v) + &
              vertical_mask(v, i) * mask(m, i)
          end do
        end do

        if (regions % function_twod == MUL_FUNC) then
          call mpas_pool_get_array(block % allFields, &
            regions % weights_twod, twod_weights, 1)

          do i = 1, solve
            do v = 1, levels
              weight_sum_oned(v) = weight_sum_oned(v) + &
                vertical_mask(v, i) * mask(m, i) * twod_weights(v, i)
            end do
          end do
        end if
      end if

      block => block % next
    end do

    ! reduce the weights and sums
    call mpas_dmpar_sum_int(domain % dminfo, count_sum_zerod, int_copy)
    count_sum_zerod = int_copy

    if (regions % function_oned == MUL_FUNC) then
      call mpas_dmpar_sum_real(domain % dminfo, weight_sum_zerod, real_copy)
      weight_sum_zerod = real_copy
    end if

    if (active_vertical) then
      call mpas_dmpar_sum_int_array(domain % dminfo, size(count_sum_oned), &
        count_sum_oned, int_vert_copy)
      count_sum_oned = int_vert_copy

      if (regions % function_twod == MUL_FUNC) then
        call mpas_dmpar_sum_real_array( &
          domain % dminfo, size(weight_sum_oned), &
          weight_sum_oned, real_vert_copy)
        weight_sum_oned = real_vert_copy
      end if
    end if

  end do

  ! free up memory
  if (active_vertical) then
    deallocate(int_vert_copy)
    deallocate(real_vert_copy)
  end if

  ! do all region reductions for each variable
  do v = 1, regions % number_of_variables
    call typed_operate(domain % dminfo, domain % blocklist, &
      regions, regions % variables(v))
  end do

  ! clean up the instance memory
  do v = 1, regions % number_of_variables
    deallocate(regions % variables(v) % output_names)
  end do
  deallocate(regions % variables)
  deallocate(regions % count_zerod_names)
  deallocate(regions % weight_zerod_names)
  deallocate(regions % count_oned_names)
  deallocate(regions % weight_oned_names)
end subroutine ocn_compute_regional_stats!}}}



!***********************************************************************
! routine ocn_restart_regional_stats
!
!> \brief Save restart for MPAS-Ocean analysis member
!> \author  Jon Woodring
!> \date    December 17, 2015
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!-----------------------------------------------------------------------
subroutine ocn_restart_regional_stats(domain, instance, err)!{{{
  ! input variables
  character (len=StrKIND), intent(in) :: instance

  ! input/output variables
  type (domain_type), intent(inout) :: domain

  ! output variables
  integer, intent(out) :: err !< Output: error flag

  ! local variables

  ! start procedure
  err = 0

end subroutine ocn_restart_regional_stats!}}}



!***********************************************************************
! routine ocn_finalize_regional_stats
!
!> \brief Finalize MPAS-Ocean analysis member
!> \author  Jon Woodring
!> \date    December 17, 2015
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!-----------------------------------------------------------------------
subroutine ocn_finalize_regional_stats(domain, instance, err)!{{{
  ! input variables
  character (len=StrKIND), intent(in) :: instance

  ! input/output variables
  type (domain_type), intent(inout) :: domain

  ! output variables
  integer, intent(out) :: err !< Output: error flag

  ! local variables

  ! start procedure
  err = 0

end subroutine ocn_finalize_regional_stats!}}}

!
! local subroutines
!

!***********************************************************************
! routine add_new_string
!
!> \brief Allocate an string in the MPAS framework for this AM
!> \author  Jon Woodring
!> \date    December 17, 2015
!> \details
!> Allocate a new integer in the AM pool and return a pointer to it.
!-----------------------------------------------------------------------
subroutine add_new_string(all_fields, inpool, outpool, field_name, target_ptr)
  ! input variables
  character (len=StrKIND) :: field_name

  ! input/output variables
  type (mpas_pool_type), pointer, intent(inout) :: all_fields, inpool, outpool

  ! output variables
  character (len=StrKIND), pointer, optional :: target_ptr

  ! local variables
  type (field0DChar), pointer :: srcString, dstString

  call mpas_pool_get_field(inpool, ONE_STRING_MEMORY, srcString, 1)
  call mpas_duplicate_field(srcString, dstString)
  dstString % fieldName = field_name
  call mpas_pool_add_field(outpool, dstString % fieldName, dstString)
  call mpas_pool_add_field(all_fields, dstString % fieldName, dstString)
  if (present(target_ptr)) then
    call mpas_pool_get_array(outpool, dstString % fieldName, target_ptr, 1)
  end if
end subroutine add_new_string



!***********************************************************************
! routine add_new_integer
!
!> \brief Allocate an integer in the MPAS framework for this AM
!> \author  Jon Woodring
!> \date    December 17, 2015
!> \details
!> Allocate a new integer in the AM pool and return a pointer to it.
!-----------------------------------------------------------------------
subroutine add_new_integer(all_fields, inpool, outpool, field_name, target_ptr)
  ! input variables
  character (len=StrKIND) :: field_name

  ! input/output variables
  type (mpas_pool_type), pointer, intent(inout) :: all_fields, inpool, outpool

  ! output variables
  integer, pointer, optional :: target_ptr

  ! local variables
  type (field0DInteger), pointer :: srcInteger, dstInteger

  call mpas_pool_get_field(inpool, ONE_INTEGER_MEMORY, srcInteger, 1)
  call mpas_duplicate_field(srcInteger, dstInteger)
  dstInteger % fieldName = field_name
  call mpas_pool_add_field(outpool, dstInteger % fieldName, dstInteger)
  call mpas_pool_add_field(all_fields, dstInteger % fieldName, dstInteger)
  if (present(target_ptr)) then
    call mpas_pool_get_array(outpool, dstInteger % fieldName, target_ptr, 1)
  end if
end subroutine add_new_integer



!***********************************************************************
! routine add_new_real
!
!> \brief Allocate an real in the MPAS framework for this AM
!> \author  Jon Woodring
!> \date    December 17, 2015
!> \details
!> Allocate a new real in the AM pool and return a pointer to it.
!-----------------------------------------------------------------------
subroutine add_new_real(all_fields, inpool, outpool, field_name, target_ptr)
  ! input variables
  character (len=StrKIND) :: field_name

  ! input/output variables
  type (mpas_pool_type), pointer, intent(inout) :: all_fields, inpool, outpool

  ! output variables
  real (kind=RKIND), pointer, optional :: target_ptr

  ! local variables
  type (field0DReal), pointer :: srcReal, dstReal

  call mpas_pool_get_field(inpool, ONE_REAL_MEMORY, srcReal, 1)
  call mpas_duplicate_field(srcReal, dstReal)
  dstReal % fieldName = field_name
  call mpas_pool_add_field(outpool, dstReal % fieldName, dstReal)
  call mpas_pool_add_field(all_fields, dstReal % fieldName, dstReal)
  if (present(target_ptr)) then
    call mpas_pool_get_array(outpool, dstReal % fieldName, target_ptr, 1)
  end if
end subroutine add_new_real



!***********************************************************************
! routine add_new_real_1d
!
!> \brief Allocate an 1d real in the MPAS framework for this AM
!> \author  Jon Woodring
!> \date    December 17, 2015
!> \details
!> Allocate a new 1d real from a 2d integer (for creating vertical
!> weight arrays from the vertical mask, i.e., drops the element dimension)
!> has_vertical is integer so that the result can be stored in the framework.
!-----------------------------------------------------------------------
subroutine add_new_real_1d(all_fields, inpool, outpool, &
    inname, outname, elem_name, &
    has_vertical, vertical_dim)!{{{
#include "regional_stats_inc/regional_field_nd_1.inc"

  type (field2DInteger), pointer :: src
  type (field1DReal), pointer :: dst
  integer, dimension(2) :: src_dims
  type (mpas_pool_field_info_type) :: info

#include "regional_stats_inc/regional_field_nd_2.inc"
  allocate(dst % array(src_dims(1)))
#include "regional_stats_inc/regional_field_nd_3.inc"

  dst % hasTimeDimension = .true.

  call mpas_pool_get_field_info(outpool, outname, info)
  info % nTimeLevels = 1
end subroutine add_new_real_1d!}}}



!***********************************************************************
! routine add_new_integer_1d
!
!> \brief Allocate an 1d integer in the MPAS framework for this AM
!> \author  Jon Woodring
!> \date    December 17, 2015
!> \details
!> Allocate a new 1d integer from a 2d integer (for creating vertical
!> count arrays from the vertical mask, i.e., drops the element dimension)
!> has_vertical is integer so that the result can be stored in the framework.
!-----------------------------------------------------------------------
subroutine add_new_integer_1d(all_fields, inpool, outpool, &
    inname, outname, elem_name, &
    has_vertical, vertical_dim)!{{{
#include "regional_stats_inc/regional_field_nd_1.inc"

  type (field2DInteger), pointer :: src
  type (field1DInteger), pointer :: dst
  integer, dimension(2) :: src_dims
  type (mpas_pool_field_info_type) :: info

#include "regional_stats_inc/regional_field_nd_2.inc"
  allocate(dst % array(src_dims(1)))
#include "regional_stats_inc/regional_field_nd_3.inc"

  dst % hasTimeDimension = .true.

  call mpas_pool_get_field_info(outpool, outname, info)
  info % nTimeLevels = 1
end subroutine add_new_integer_1d!}}}



!***********************************************************************
! function output_naming
!
!> \brief Given an input name, create a corresponding output name
!> \author  Jon Woodring
!> \date    December 17, 2015
!> \details
!> Code to create consistent output names from input names.
!-----------------------------------------------------------------------
character (len=StrKIND) function output_naming &
(op_name, input_name, region_identifier, instance)
  character (len=StrKIND), intent(in) :: op_name, &
    input_name, region_identifier, instance

  output_naming = trim(instance) // &
    trim(region_identifier) // '_' // trim(op_name) // '_' // &
    trim(input_name)
end function output_naming



!***********************************************************************
! function operator_naming
!
!> \brief Given an operator enum, create a corresponding operator name
!> \author  Jon Woodring
!> \date    December 17, 2015
!> \details
!> Code to create consistent operator names from operator enums.
!-----------------------------------------------------------------------
character (len=StrKIND) function operator_naming(enum)
  integer, intent(in) :: enum

  if (enum == AVG_OP) then
    operator_naming = AVG_TOKEN
  else if (enum == MIN_OP) then
    operator_naming = MIN_TOKEN
  else if (enum == MAX_OP) then
    operator_naming = MAX_TOKEN
  else if (enum == SUM_OP) then
    operator_naming = SUM_TOKEN
  else if (enum == SOS_OP) then
    operator_naming = SOS_TOKEN
  else
    call mpas_dmpar_global_abort( &
      trim(CURRENT_CORE_NAME) // ' ERROR: The impossible happened. ' // &
      'Tried creating an operator of unknown kind in regional stats AM.')
  end if
end function operator_naming



!***********************************************************************
! function element_naming
!
!> \brief Given an element enum, create a corresponding element name
!> \author  Jon Woodring
!> \date    December 17, 2015
!> \details
!> Code to create consistent operator names from operator enums.
!-----------------------------------------------------------------------
character (len=StrKIND) function element_naming(enum)
  integer, intent(in) :: enum

  if (enum == CELL_REGION) then
    element_naming = CELL_TOKEN
  else if (enum == VERTEX_REGION) then
    element_naming = VERTEX_TOKEN
  else
    call mpas_dmpar_global_abort( &
      trim(CURRENT_CORE_NAME) // ' ERROR: The impossible happened. ' // &
      'Tried creating an element of unknown kind in regional stats AM.')
  end if
end function element_naming



!***********************************************************************
! function dimension_naming
!
!> \brief Given a region type, return its canonical dimension string name
!> \author  Jon Woodring
!> \date    December 17, 2015
!> \details
!> Returns the canonical MPAS dimension name given a region element dimension.
!-----------------------------------------------------------------------
character (len=StrKIND) function dimension_naming(elem_type)
  integer, intent(in) :: elem_type

  if (elem_type == CELL_REGION) then
    dimension_naming = CELL_DIM
  else if (elem_type == VERTEX_REGION) then
    dimension_naming = VERTEX_DIM
  else
    call mpas_dmpar_global_abort( &
      trim(CURRENT_CORE_NAME) // ' ERROR: The impossible happened. ' // &
      'Tried creating a element dimension name of unknown kind in ' // &
      'regional stats AM.')
  end if
end function dimension_naming



!***********************************************************************
! function check_real_element_dim
!
!> \brief Return true if the last dimension matches the element and is real.
!> \author  Jon Woodring
!> \date    December 17, 2015
!> \details
!> Given an element dimension name and field info, return true
!> if the last dimension of the field is the element and the field is
!> real, otherwise return false.
!-----------------------------------------------------------------------
logical function check_real_element_dim(all_fields, field_name, elem_name)
  type (mpas_pool_type), pointer, intent(in) :: all_fields
  character (len=StrKIND), intent(in) :: field_name, elem_name

  logical :: has_elem
  type (mpas_pool_field_info_type) :: info
  type (field1DReal), pointer :: r1
  type (field2DReal), pointer :: r2
  type (field3DReal), pointer :: r3
  type (field4DReal), pointer :: r4
  type (field5DReal), pointer :: r5

  call mpas_pool_get_field_info(all_fields, field_name, info)

  if (info % fieldType /= MPAS_POOL_REAL) then
    write(stdErrUnit,*) &
      trim(CURRENT_CORE_NAME) // ' WARNING: field "' // &
      trim(field_name) // '" listed in the ' // &
      'output stream, for regional stats analysis member, ' // &
      'is not real. Regional stats will not be applied to ' // &
      'this field.'

    check_real_element_dim = .false.
  else if (info % nDims < 1) then
    write(stdErrUnit,*) &
      trim(CURRENT_CORE_NAME) // ' WARNING: field "' // &
      trim(field_name) // '" listed in the ' // &
      'output stream, for regional stats analysis member, ' // &
      'is 0D. Regional stats will not be applied to ' // &
      'this field.'

    check_real_element_dim = .false.
  else
    if (info % nDims == 1) then
      call mpas_pool_get_field(all_fields, field_name, r1, 1)
      has_elem = check_element_dim(r1 % dimNames, elem_name)
    else if (info % nDims == 2) then
      call mpas_pool_get_field(all_fields, field_name, r2, 1)
      has_elem = check_element_dim(r2 % dimNames, elem_name)
    else if (info % nDims == 3) then
      call mpas_pool_get_field(all_fields, field_name, r3, 1)
      has_elem = check_element_dim(r3 % dimNames, elem_name)
    else if (info % nDims == 4) then
      call mpas_pool_get_field(all_fields, field_name, r4, 1)
      has_elem = check_element_dim(r4 % dimNames, elem_name)
    else
      call mpas_pool_get_field(all_fields, field_name, r5, 1)
      has_elem = check_element_dim(r5 % dimNames, elem_name)
    end if

    if (.not. has_elem) then
      write(stdErrUnit,*) &
        trim(CURRENT_CORE_NAME) // ' WARNING: field "' // &
        trim(field_name) // '" listed in the output stream, ' // &
        'for regional stats analysis member, does not have ' // &
        trim(elem_name) // ' as its last dimension. Regional stats ' // &
        'will not be applied to this field.'
    end if

    check_real_element_dim = has_elem
  end if

end function check_real_element_dim



!***********************************************************************
! function check_element_dim
!
!> \brief Return true if the last dimension matches the element name.
!> \author  Jon Woodring
!> \date    December 17, 2015
!> \details
!> Given an element dimension name and list of dim names, return true
!> if the last dimension is the element, otherwise return false.
!-----------------------------------------------------------------------
logical function check_element_dim(dim_names, elem_name)
  character (len=StrKIND), dimension(:), intent(in) :: dim_names
  character (len=StrKIND), intent(in) :: elem_name

  integer :: last

  last = size(dim_names)
  if (last > 0) then
    check_element_dim = (trim(dim_names(last)) == trim(elem_name))
  else
    check_element_dim = .false.
  end if
end function check_element_dim



!***********************************************************************
! function check_vertical_dim
!
!> \brief Return not 0 if the last dimension matches the vertical name.
!> \author  Jon Woodring
!> \date    December 17, 2015
!> \details
!> Given an vertical dimension name and list of dim names, return not 0
!> if the last dimension is the vertical, otherwise return 0.
!> It uses integers so the result can be stored in the framework.
!-----------------------------------------------------------------------
integer function check_vertical_dim(dim_names, vert_name)
  character (len=StrKIND), dimension(:), intent(in) :: dim_names
  character (len=StrKIND), intent(in) :: vert_name

  integer :: last

  last = size(dim_names)
  if (last > 1) then
    if (trim(dim_names(last - 1)) == trim(vert_name)) then
      check_vertical_dim = 1
    else
      check_vertical_dim = 0
    end if
  else
    check_vertical_dim = 0
  end if
end function check_vertical_dim



!***********************************************************************
! function mask_naming
!
!> \brief Consistent naming for mask field from the stream.
!> \author  Jon Woodring
!> \date    December 17, 2015
!> \details Creates the name of the mask field based on the element type.
!-----------------------------------------------------------------------
character (len=StrKIND) function mask_naming(elem_type)
  integer, intent(in) :: elem_type

  if (elem_type == CELL_REGION) then
    mask_naming = &
      trim(MASK_DATA_PREFIX) // trim(CELL_CAMEL) // trim(MASK_DATA_SUFFIX)
  else if (elem_type == VERTEX_REGION) then
    mask_naming = &
      trim(MASK_DATA_PREFIX) // trim(VERTEX_CAMEL) // trim(MASK_DATA_SUFFIX)
  else
    call mpas_dmpar_global_abort( &
      trim(CURRENT_CORE_NAME) // ' ERROR: The impossible happened. ' // &
      'Tried to create a name for a mask field of unknown type.')
  end if
end function mask_naming


!***********************************************************************
! function fix_region_name
!
!> \brief Replaces non identifier characters with _ (underscore)
!> \author  Jon Woodring
!> \date    December 17, 2015
!> \details Given a character string, it will replace all of the Fortran
!  non-identifier characters (anything but [a-zA-Z0-9_]) with underscore.
!  This is primarily for renaming regions from the mask file,
!  because NetCDF doesn't like non-identifier characters in the
!  start position and Fortran identifiers are [a-zA-Z][a-zA-Z0-9_]+.
!  It will also prepend 'reg' to the string to make sure it starts
!  with [a-zA-Z].
!-----------------------------------------------------------------------
character (len=StrKIND) function fix_region_name(string)
  character (len=StrKIND), intent(in) :: string

  character (len=StrKIND) :: copy
  integer :: tl, i, c

  copy = string
  tl = len_trim(copy)

  if (tl < 1) then
    call mpas_dmpar_global_abort( &
      trim(CURRENT_CORE_NAME) // ' ERROR: a region name, from the  ' // &
      'region mask input for the regional stats analysis member, the ' // &
      'is an empty string.')
  end if

  tl = tl + 3
  copy = 'reg' // copy

  do i = 4, tl
    c = iachar(copy(i:i))
    if ((c < 48) .or. ((c > 57) .and. (c < 65)) .or. (c > 122)) then
      copy(i:i) = '_'
    end if
  end do

  fix_region_name = copy
end function fix_region_name


!***********************************************************************
! routine debug_state
!
!> \brief Print all of the state for this instance.
!> \author  Jon Woodring
!> \date    December 17, 2015
!> \details
!> Given the current state in regions, print it all for debugging.
!-----------------------------------------------------------------------
subroutine debug_state(regions)
  type (regional_type), intent(in) :: regions
  integer :: v, b, last, m

  write(*,*) 'num vars: ', regions % number_of_variables

  write (*,*) 'op: ', regions % operation
  write (*,*) 'elem: ', regions % region_element
  write (*,*) '1d func: ', regions % function_oned
  write (*,*) '1d weight: ', trim(regions % weights_oned)
  write (*,*) '2d func: ', regions % function_twod
  write (*,*) '2d weight: ', trim(regions % weights_twod)

  write (*,*) 'mask: ', trim(regions % masking_field)
  write (*,*) 'vertical mask: ', trim(regions % vertical_mask)

  last = regions % num_regions_per(regions % group_index)
  write(*,*) 'selected group: ', regions % group_index
  write(*,*) 'num regions in: ', last
  do b = 1, last
    write(*,*) 'count 1d var: ', trim(regions % count_zerod_names(b))
    write(*,*) 'weight 1d var: ', trim(regions % weight_zerod_names(b))
    write(*,*) 'count 2d var: ', trim(regions % count_oned_names(b))
    write(*,*) 'weight 2d var: ', trim(regions % weight_oned_names(b))
  end do
  do v = 1, regions % number_of_variables
    write(*,*) 'var: ', trim(regions % variables(v) % input_name)
    write(*,*) 'has vertical: ', regions % variables(v) % has_vertical
    do b = 1, last
      write(*,*) 'index: ', b
      m = regions % groups(b, regions % group_index)
      write(*,*) 'region index: ', m
      write(*,*) 'out var: ', trim(regions % variables(v) % output_names(b))
    end do
  end do
end subroutine debug_state



!***********************************************************************
! routine get_state
!
!> \brief Get all of the state for this instance.
!> \author  Jon Woodring
!> \date    December 17, 2015
!> \details
!> This will allocate and fetch all of the state necessary for this
!> instance that is being run (prevalidated data from start_state).
!-----------------------------------------------------------------------
subroutine get_state(domain, instance, regions)
  ! input variables
  character (len=StrKIND), intent(in) :: instance

  ! input/output variables
  type (domain_type), intent(inout) :: domain

  ! output variables
  type (regional_type), intent(out) :: regions

  ! local variables
  integer :: v, b, last, m
  character (len=StrKIND) :: namelist_prefix, storage_prefix, field_name, &
    var_identifier, op_name
  type (mpas_pool_type), pointer :: maskPool, amPool
  character (len=StrKIND), dimension(:), pointer :: names
  logical :: active_vertical

  ! start procedure
  namelist_prefix = trim(CONFIG_PREFIX) // trim(instance)
  storage_prefix = trim(FRAMEWORK_PREFIX) // trim(instance)

  call mpas_pool_get_subpool(domain % blocklist % structs, &
    REGIONAL_STATS_POOL, amPool)

  !
  ! get ones we stored in framework for this instance
  !

  ! number_of_variables
  field_name = trim(storage_prefix) // trim(NUMBER_OF_VARIABLES_SUFFIX)
  call mpas_pool_get_array(amPool, field_name, regions % number_of_variables, 1)

  ! operation
  field_name = trim(storage_prefix) // trim(OPERATION_SUFFIX)
  call mpas_pool_get_array(amPool, field_name, regions % operation, 1)

  op_name = operator_naming(regions % operation)

  ! mask type
  field_name = trim(storage_prefix) // trim(REGION_TYPE_SUFFIX)
  call mpas_pool_get_array(amPool, field_name, regions % region_element, 1)

  ! mask field & pool
  regions % masking_field = mask_naming(regions % region_element)
  call mpas_pool_get_subpool(domain % blocklist % structs, &
    MASK_POOL_NAME, maskPool)

  ! weighting function
  field_name = trim(storage_prefix) // trim(FUNCTION_ONED_SUFFIX)
  call mpas_pool_get_array(amPool, field_name, regions % function_oned, 1)

  field_name = trim(storage_prefix) // trim(FUNCTION_TWOD_SUFFIX)
  call mpas_pool_get_array(amPool, field_name, regions % function_twod, 1)

  ! get the input names and has vertical for variables
  allocate(regions % variables(regions % number_of_variables))
  do v = 1, regions % number_of_variables
    ! identifier
    write(var_identifier, '(I0)') v

    ! input name
    field_name = trim(storage_prefix) // '_' // trim(var_identifier) // &
      trim(INPUT_NAME_SUFFIX)
    call mpas_pool_get_array(amPool, field_name, &
      regions % variables(v) % input_name, 1)

    ! has_vertical
    field_name = trim(storage_prefix) // '_' // trim(var_identifier) // &
      trim(HAS_VERTICAL_SUFFIX)
    call mpas_pool_get_array(amPool, field_name, &
      regions % variables(v) % has_vertical, 1)
  end do

  ! group index
  field_name = trim(storage_prefix) // trim(REGION_GROUP_SUFFIX)
  call mpas_pool_get_array(amPool, field_name, regions % group_index, 1)

  !
  ! get ones that already exist from namelist/stream
  !

  ! weighting fields
  field_name = trim(namelist_prefix) // trim(WEIGHTS_ONED_SUFFIX)
  call mpas_pool_get_config(domain % configs, field_name, &
    regions % weights_oned)

  field_name = trim(namelist_prefix) // trim(WEIGHTS_TWOD_SUFFIX)
  call mpas_pool_get_config(domain % configs, field_name, &
    regions % weights_twod)

  ! groups
  call mpas_pool_get_array(maskPool, GROUP_DATA_NAME, regions % groups, 1)

  ! num_regions_per
  call mpas_pool_get_array(maskPool, REGIONS_PER_NAME, &
    regions % num_regions_per, 1)

  ! vertical mask
  field_name = trim(namelist_prefix) // trim(VERTICAL_MASK_SUFFIX)
  call mpas_pool_get_config(domain % configs, field_name, &
    regions % vertical_mask)

  field_name = trim(namelist_prefix) // trim(VERTICAL_DIM_SUFFIX)
  call mpas_pool_get_config(domain % configs, field_name, &
    regions % vertical_dim)

  active_vertical = (trim(regions % vertical_mask) /= trim(NONE_TOKEN))

  !
  ! generate output names
  !

  ! fetch region names and name via regions
  call mpas_pool_get_array(maskPool, REGION_NAMES_NAME, names)

  last = regions % num_regions_per(regions % group_index)
  ! allocate count and weight names
  allocate(regions % count_zerod_names(last))
  allocate(regions % weight_zerod_names(last))
  allocate(regions % count_oned_names(last))
  allocate(regions % weight_oned_names(last))
  do b = 1, last
    m = regions % groups(b, regions % group_index)

    ! 0d count & weight name
    regions % count_zerod_names(b) = output_naming(COUNT_TOKEN, &
      ZEROD_TOKEN, fix_region_name(names(m)), instance)

    if (regions % function_oned == MUL_FUNC) then
      regions % weight_zerod_names(b) = &
        output_naming(WEIGHT_TOKEN, &
        ZEROD_TOKEN, fix_region_name(names(m)), instance)
    end if

    ! 1d count & weight name
    if (active_vertical) then
      regions % count_oned_names(b) = &
        output_naming(COUNT_TOKEN, &
        ONED_TOKEN, fix_region_name(names(m)), instance)

      if (regions % function_twod == MUL_FUNC) then
        regions % weight_oned_names(b) = &
          output_naming(WEIGHT_TOKEN, &
          ONED_TOKEN, fix_region_name(names(m)), instance)
      end if
    end if

  end do

  do v = 1, regions % number_of_variables
    ! allocate output names
    allocate(regions % variables(v) % output_names(last))

    do b = 1, last
      m = regions % groups(b, regions % group_index)

      ! region name
      regions % variables(v) % output_names(b) = output_naming(op_name, &
        regions % variables(v) % input_name, &
        fix_region_name(names(m)), instance)
    end do
  end do

end subroutine get_state


!***********************************************************************
! routine start_state
!
!> \brief Begin the initialization of this analysis member
!> \author  Jon Woodring
!> \date    December 17, 2015
!> \details
!> All the necessary details to initialize this analysis member
!> instance.
!-----------------------------------------------------------------------
subroutine start_state(domain, instance, err)
  ! input variables
  character (len=StrKIND), intent(in) :: instance

  ! input/output variables
  type (domain_type), intent(inout) :: domain

  ! output variables
  integer, intent(out) :: err !< Output: error flag

  ! local variables
  type (regional_type) :: regions
  character (len=StrKIND), pointer :: config_results, output_stream_name, &
    config_results_2
  character (len=StrKIND) :: namelist_prefix, storage_prefix, field_name, &
    elem_name, op_name, var_identifier
  integer :: v, b, last, m
  integer, pointer :: number_of_regions, number_of_groups, max_regions_per
  type (mpas_pool_type), pointer :: maskPool, amPool
  type (field1DReal), pointer :: oned_field
  type (field2DReal), pointer :: twod_field
  type (field2DInteger), pointer :: mask
  type (mpas_pool_field_info_type) :: info
  character (len=StrKIND), dimension(:), pointer :: names
  logical :: active_vertical
  logical, dimension(:), allocatable :: valid_input

  ! start procedure
  err = 0

  namelist_prefix = trim(CONFIG_PREFIX) // trim(instance)
  storage_prefix = trim(FRAMEWORK_PREFIX) // trim(instance)

  call mpas_pool_get_subpool(domain % blocklist % structs, &
    REGIONAL_STATS_POOL, amPool)

  !
  ! get dimensions
  !

  call mpas_pool_get_dimension(domain % blocklist % dimensions, &
    NUM_REGIONS_SUFFIX, number_of_regions)
  call mpas_pool_get_dimension(domain % blocklist % dimensions, &
    NUM_GROUPS_SUFFIX, number_of_groups)
  call mpas_pool_get_dimension(domain % blocklist % dimensions, &
    MAX_REGIONS_SUFFIX, max_regions_per)

  if (number_of_regions < 1) then
    call mpas_dmpar_global_abort( &
      trim(CURRENT_CORE_NAME) // ' ERROR: number of regions dimension is ' // &
      'less than 1 for the regional stats analysis member.')
  end if

  if (number_of_groups < 1) then
    call mpas_dmpar_global_abort( &
      trim(CURRENT_CORE_NAME) // ' ERROR: number of region groups ' // &
      'dimension is less than 1 for the regional stats analysis member.')
  end if

  if (max_regions_per < 1) then
    call mpas_dmpar_global_abort( &
      trim(CURRENT_CORE_NAME) // ' ERROR: maximum number of regions per ' // &
      'group dimension is less than 1 for the regional stats analysis member.')
  end if

  !
  ! allocate some framework memory for instance state, and assign values
  !

  ! number of variables done in modify stream

  ! operation
  field_name = trim(storage_prefix) // trim(OPERATION_SUFFIX)
  call add_new_integer(domain % blocklist % allFields, amPool, amPool, &
    field_name, regions % operation)

  ! get our operation
  field_name = trim(namelist_prefix) // trim(OPERATION_SUFFIX)
  call mpas_pool_get_config(domain % configs, field_name, config_results)
  if (config_results == AVG_TOKEN) then
    regions % operation = AVG_OP
  else if (config_results == MIN_TOKEN) then
    regions % operation = MIN_OP
  else if (config_results == MAX_TOKEN) then
    regions % operation = MAX_OP
  else if (config_results == SUM_TOKEN) then
    regions % operation = SUM_OP
  else if (config_results == SOS_TOKEN) then
    regions % operation = SOS_OP
  else
    ! error if unknown operation
    call mpas_dmpar_global_abort( &
      trim(CURRENT_CORE_NAME) // ' ERROR: unknown operation in ' // &
      'regional stats analysis member configuration.')
  end if

  op_name = operator_naming(regions % operation)

  ! element type
  field_name = trim(storage_prefix) // trim(REGION_TYPE_SUFFIX)
  call add_new_integer(domain % blocklist % allFields, amPool, amPool, &
    field_name, regions % region_element)

  ! get the element type
  field_name = trim(namelist_prefix) // trim(REGION_TYPE_SUFFIX)
  call mpas_pool_get_config(domain % configs, field_name, config_results)

  if (config_results == CELL_TOKEN) then
    regions % region_element = CELL_REGION
  else if (config_results == VERTEX_TOKEN) then
    regions % region_element = VERTEX_REGION
  else
    ! error if unknown element type
    call mpas_dmpar_global_abort( &
      trim(CURRENT_CORE_NAME) // ' ERROR: unknown region element type in ' // &
      'regional stats analysis member configuration.')
  end if

  elem_name = dimension_naming(regions % region_element)

  ! get vertical dim
  field_name = trim(namelist_prefix) // trim(VERTICAL_DIM_SUFFIX)
  call mpas_pool_get_config(domain % configs, field_name, &
    regions % vertical_dim)

  active_vertical = (trim(regions % vertical_dim) /= trim(NONE_TOKEN))

  ! mask field & pool
  regions % masking_field = mask_naming(regions % region_element)
  call mpas_pool_get_subpool(domain % blocklist % structs, &
    MASK_POOL_NAME, maskPool)

  ! validate the mask field
  !
  ! shouldn't actually need this because the naming convention, but
  ! can't hurt to check, just in case
  if (regions % masking_field /= NONE_TOKEN) then
    call mpas_pool_get_field_info(maskPool, regions % masking_field, info)

    if (info % fieldType /= MPAS_POOL_INTEGER) then
      call mpas_dmpar_global_abort( &
        trim(CURRENT_CORE_NAME) // ' ERROR: the mask ' // &
        trim(regions % masking_field) // &
        ' for the regional stats AM is not a integer field.')
    end if

    if (info % nDims /= 2) then
      call mpas_dmpar_global_abort( &
        trim(CURRENT_CORE_NAME) // ' ERROR: the mask ' // &
        trim(regions % masking_field) // &
        ' for the regional stats AM is not a 2D field.')
    end if

    call mpas_pool_get_field(maskPool, regions % masking_field, mask)

    if (.not. check_element_dim(mask % dimNames, elem_name)) then
      call mpas_dmpar_global_abort( &
        trim(CURRENT_CORE_NAME) // ' ERROR: the mask ' // &
        trim(regions % masking_field) // &
        ' for the regional stats AM needs to have ' // trim(elem_name) // &
        ' as its last dimension.')
    end if
  end if

  ! 1d weighting function
  field_name = trim(storage_prefix) // trim(FUNCTION_ONED_SUFFIX)
  call add_new_integer(domain % blocklist % allFields, amPool, amPool, &
    field_name, regions % function_oned)

  field_name = trim(namelist_prefix) // trim(FUNCTION_ONED_SUFFIX)
  call mpas_pool_get_config(domain % configs, field_name, config_results)

  if (config_results == ID_TOKEN) then
    regions % function_oned = ID_FUNC
  else if (config_results == MUL_TOKEN) then
    regions % function_oned = MUL_FUNC
  else
    ! error if unknown operation
    call mpas_dmpar_global_abort( &
      trim(CURRENT_CORE_NAME) // &
      ' ERROR: unknown 1d weighting function in ' // &
      'regional stats analysis member configuration.')
  end if

  if (regions % function_oned == MUL_FUNC .and. &
      (regions % operation == MIN_OP .or. regions % operation == MAX_OP)) then
    call mpas_dmpar_global_abort( &
      trim(CURRENT_CORE_NAME) // ' ERROR: 1d weighting function in ' // &
      'regional stats analysis member can only "mul" if the' // &
      'if the operation is "avg", "sum", or "sos"')
  end if

  ! 2d weighting function
  field_name = trim(storage_prefix) // trim(FUNCTION_TWOD_SUFFIX)
  call add_new_integer(domain % blocklist % allFields, amPool, amPool, &
    field_name, regions % function_twod)

  field_name = trim(namelist_prefix) // trim(FUNCTION_TWOD_SUFFIX)
  call mpas_pool_get_config(domain % configs, field_name, config_results)

  if (config_results == ID_TOKEN) then
    regions % function_twod = ID_FUNC
  else if (config_results == MUL_TOKEN) then
    regions % function_twod = MUL_FUNC
  else
    ! error if unknown operation
    call mpas_dmpar_global_abort( &
      trim(CURRENT_CORE_NAME) // &
      ' ERROR: unknown 2d weighting function in ' // &
      'regional stats analysis member configuration.')
  end if

  if (regions % function_twod == MUL_FUNC .and. &
      (regions % operation == MIN_OP .or. regions % operation == MAX_OP)) then
    call mpas_dmpar_global_abort( &
      trim(CURRENT_CORE_NAME) // ' ERROR: 2d weighting function in ' // &
      'regional stats analysis member can only "mul" if the' // &
      'if the operation is "avg", "sum", or "sos"')
  end if

  ! group index
  field_name = trim(storage_prefix) // trim(REGION_GROUP_SUFFIX)
  call add_new_integer(domain % blocklist % allFields, amPool, amPool, &
    field_name, regions % group_index)

  field_name = trim(namelist_prefix) // trim(REGION_GROUP_SUFFIX)
  call mpas_pool_get_config(domain % configs, field_name, &
    config_results)

  call mpas_pool_get_array(maskPool, GROUP_NAMES_NAME, names)

  ! error if we can't find it
  regions % group_index = 0
  do v = 1, number_of_groups
    if (trim(names(v)) == trim(config_results)) then
      regions % group_index = v
    end if
  end do

  if (regions % group_index == 0) then
    call mpas_dmpar_global_abort( &
      trim(CURRENT_CORE_NAME) // ' ERROR: unable to find the region group ' // &
      'named "' // trim(config_results) // '" in the region mask input stream')
  end if

  !
  ! verify the other data is OK before trying to use it later
  !

  ! validate input stream is not none
  field_name = trim(namelist_prefix) // trim(INPUT_STREAM_SUFFIX)
  call mpas_pool_get_config(domain % configs, field_name, config_results)

  if (config_results == NONE_TOKEN) then
    call mpas_dmpar_global_abort( &
      trim(CURRENT_CORE_NAME) // ' ERROR: input stream for regional ' // &
      'stats AM cannot be "none". It needs to point to the region/mask stream.')
  end if

  ! validate the restart stream is not none
  field_name = trim(namelist_prefix) // trim(RESTART_STREAM_SUFFIX)
  call mpas_pool_get_config(domain % configs, field_name, config_results_2)

  if (config_results == NONE_TOKEN) then
    call mpas_dmpar_global_abort( &
      trim(CURRENT_CORE_NAME) // ' ERROR: restart stream for regional ' // &
      'stats AM cannot be "none". It needs to point to the region/mask stream.')
  end if

  ! make sure input and restart are the same
  if (config_results /= config_results_2) then
    call mpas_dmpar_global_abort( &
      trim(CURRENT_CORE_NAME) // ' ERROR: the input and restart stream ' // &
      'for regional stats AM need to be the same to ensure that it ' // &
      'has the same behavior on start and restart.')
  end if

  ! validate vertical mask
  field_name = trim(namelist_prefix) // trim(VERTICAL_MASK_SUFFIX)
  call mpas_pool_get_config(domain % configs, field_name, &
    regions % vertical_mask)

  active_vertical = .false.
  if (regions % vertical_mask /= NONE_TOKEN) then
    call mpas_pool_get_field_info(domain % blocklist % allFields, &
      regions % vertical_mask, info)

    if (info % fieldType /= MPAS_POOL_INTEGER) then
      call mpas_dmpar_global_abort( &
        trim(CURRENT_CORE_NAME) // ' ERROR: the vertical mask ' // &
        trim(regions % vertical_mask) // &
        'for the regional stats AM is not an integer field.')
    end if

    if (info % nDims /= 2) then
      call mpas_dmpar_global_abort( &
        trim(CURRENT_CORE_NAME) // ' ERROR: the vertical mask ' // &
        trim(regions % vertical_mask) // &
        ' for the regional stats AM is not a 2D field.')
    end if

    call mpas_pool_get_field(domain % blocklist % allFields, &
      regions % vertical_mask, mask)

    if (.not. check_element_dim(mask % dimNames, elem_name)) then
      call mpas_dmpar_global_abort( &
        trim(CURRENT_CORE_NAME) // ' ERROR: the vertical mask ' // &
        trim(regions % vertical_mask) // &
        ' for the regional stats AM needs to have ' // trim(elem_name) // &
        ' as its last dimension.')
    end if

    if (check_vertical_dim(mask % dimNames, regions % vertical_dim) == 0) then
      call mpas_dmpar_global_abort( &
        trim(CURRENT_CORE_NAME) // ' ERROR: the vertical mask ' // &
        trim(regions % vertical_mask) // &
        ' for the regional stats AM does not have ' // &
        trim(regions % vertical_dim) // &
        ' as its second to last dimension.')
    end if

    active_vertical = .true.
  end if

   ! validate 1d weighting field
  field_name = trim(namelist_prefix) // trim(WEIGHTS_ONED_SUFFIX)
  call mpas_pool_get_config(domain % configs, field_name, &
    regions % weights_oned)

  if (regions % function_oned == ID_FUNC) then
    if (regions % weights_oned /= NONE_TOKEN) then
      call mpas_dmpar_global_abort( &
        trim(CURRENT_CORE_NAME) // ' ERROR: 1d weighting field "' // &
        trim(regions % weights_oned) // '" is not set to "none" ' // &
        'when the 1d weighting function is set to "id" in ' // &
        'regional stats AM.')
    end if
  else
    ! weighting field & info
    if (regions % weights_oned == NONE_TOKEN) then
      call mpas_dmpar_global_abort( &
        trim(CURRENT_CORE_NAME) // ' ERROR: 1d weighting field "' // &
        trim(regions % weights_oned) // '" is set to "none" ' // &
        'when the 1d weighting function is something other than "id" in ' // &
        'regional stats AM.')
    else
      call mpas_pool_get_field_info(domain % blocklist % allFields, &
        regions % weights_oned, info)

      if (info % nDims /= 1) then
        call mpas_dmpar_global_abort( &
          trim(CURRENT_CORE_NAME) // &
          ' ERROR: the listed 1d weighting field "' // &
          trim(regions % weights_oned) // &
          '" for the regional stats AM is not actually a 1D field.')
      end if

      ! check if we can handle it
      if (info % fieldType == MPAS_POOL_REAL) then
        call mpas_pool_get_field(domain % blocklist % allFields, &
          regions % weights_oned, oned_field)

        if (.not. check_element_dim(oned_field % dimNames, elem_name)) then
          call mpas_dmpar_global_abort( &
            trim(CURRENT_CORE_NAME) // ' ERROR: 1d weighting field "' // &
            regions % weights_oned // &
            '" in regional stats AM needs to have ' // &
            elem_name // ' dimensioni as its last dimension.')
        end if
      else
        call mpas_dmpar_global_abort( &
          trim(CURRENT_CORE_NAME) // ' ERROR: 1d weighting field "' // &
          trim(regions % weights_oned) // '" listed in the ' // &
          'namelist, for regional stats analysis member ' // &
          'stream, is not real.')
      end if
    end if
  end if

  ! validate 2d weighting field
  field_name = trim(namelist_prefix) // trim(WEIGHTS_TWOD_SUFFIX)
  call mpas_pool_get_config(domain % configs, field_name, &
    regions % weights_twod)

  if (.not. active_vertical) then
    if (regions % function_twod /= ID_FUNC) then
      call mpas_dmpar_global_abort( &
        trim(CURRENT_CORE_NAME) // ' ERROR: vertical dimension is not set ' // &
        'in regional stats AM when the 2d weight function is something ' // &
        'than "id"')
    end if
  else if (regions % function_twod == ID_FUNC) then
    if (regions % weights_twod /= NONE_TOKEN) then
      call mpas_dmpar_global_abort( &
        trim(CURRENT_CORE_NAME) // ' ERROR: 2d weighting field "' // &
        trim(regions % weights_twod) // '" is not set to "none"' // &
        'when the 2d weighting function is set to "id", in ' // &
        'regional stats AM.')
    end if
  else
    ! weighting field & info
    if (regions % weights_twod == NONE_TOKEN) then
      call mpas_dmpar_global_abort( &
        trim(CURRENT_CORE_NAME) // ' ERROR: 2d weighting field "' // &
        trim(regions % weights_twod) // '" is set to "none"' // &
        'when the 2d weighting function is something other than "id" in ' // &
        'regional stats AM.')
    else
      call mpas_pool_get_field_info(domain % blocklist % allFields, &
        regions % weights_twod, info)

      if (info % nDims /= 2) then
        call mpas_dmpar_global_abort( &
          trim(CURRENT_CORE_NAME) // &
          ' ERROR: the listed 2d weighting field "' // &
          trim(regions % weights_twod) // &
          '" for the regional stats AM is not actually a 2D field.')
      end if

      ! check if we can handle it
      if (info % fieldType == MPAS_POOL_REAL) then
        call mpas_pool_get_field(domain % blocklist % allFields, &
          regions % weights_twod, twod_field)

        if (.not. check_element_dim(twod_field % dimNames, elem_name)) then
          call mpas_dmpar_global_abort( &
            trim(CURRENT_CORE_NAME) // ' ERROR: 2d weighting field "' // &
            regions % weights_twod // &
            '" in regional stats AM needs to have ' // &
            elem_name // ' as its last dimension.')
        end if

        if (check_vertical_dim(twod_field % dimNames, &
            regions % vertical_dim) == 0) then
          call mpas_dmpar_global_abort( &
            trim(CURRENT_CORE_NAME) // ' ERROR: 2d weighting field "' // &
            regions % weights_twod // &
            '" in regional stats AM needs to have ' // &
            regions % vertical_dim // &
            ' vertical dimension as its second to last dimension.')
        end if
      else
        call mpas_dmpar_global_abort( &
          trim(CURRENT_CORE_NAME) // ' ERROR: 2d weighting field "' // &
          trim(regions % weights_twod) // '" listed in the ' // &
          'namelist, for regional stats analysis member ' // &
          'stream, is not real.')
      end if
    end if
  end if

  ! groups
  call mpas_pool_get_array(maskPool, GROUP_DATA_NAME, regions % groups, 1)

  ! num_regions
  call mpas_pool_get_array(maskPool, REGIONS_PER_NAME, &
    regions % num_regions_per, 1)

  ! verify selection is OK
  last = regions % num_regions_per(regions % group_index)
  if (last > max_regions_per) then
    call mpas_dmpar_global_abort( &
      trim(CURRENT_CORE_NAME) // ' ERROR: the number of regions for ' // &
      'selected group, in the regional stats analysis member, ' // &
      'is greater than the maximum number of regions dimension.')
  end if
  if (last < 1) then
    call mpas_dmpar_global_abort( &
      trim(CURRENT_CORE_NAME) // ' ERROR: the number of regions for ' // &
      'selection group, in the regional stats analysis member, ' // &
      'is less than 1.')
  end if

  do b = 1, last
    v = regions % groups(b, regions % group_index)
    if ((v < 1) .or. (v > number_of_regions)) then
      call mpas_dmpar_global_abort( &
        trim(CURRENT_CORE_NAME) // ' ERROR: a region index found ' // &
        'in the list of regions for the selected group in the ' // &
        'regional stats analysis member is out of bounds for the size ' // &
        'of the regions found in the masks input.')
    end if
  end do


  !
  ! modify the stream to create destination output variables and remove inputs
  !

  ! get the output stream name
  field_name = trim(namelist_prefix) // trim(OUTPUT_STREAM_SUFFIX)
  call mpas_pool_get_config(domain % configs, field_name, output_stream_name)

  if (output_stream_name == NONE_TOKEN) then
    call mpas_dmpar_global_abort( &
      trim(CURRENT_CORE_NAME) // ' ERROR: output stream cannot be "none" ' // &
      'for regional stats.')
  end if

  ! number_of_variables
  field_name = trim(storage_prefix) // trim(NUMBER_OF_VARIABLES_SUFFIX)
  call add_new_integer(domain % blocklist % allFields, amPool, amPool, &
    field_name, regions % number_of_variables)

  ! count total number of fields
  call mpas_stream_mgr_begin_iteration(domain % streamManager, &
    output_stream_name, err)
  b = 0
  do while (mpas_stream_mgr_get_next_field(domain % streamManager, &
    output_stream_name, field_name))
    b = b + 1
  end do

  allocate(valid_input(b))

  ! count the number of variables and mark if valid
  call mpas_stream_mgr_begin_iteration(domain % streamManager, &
    output_stream_name, err)
  b = 1
  regions % number_of_variables = 0
  do while (mpas_stream_mgr_get_next_field(domain % streamManager, &
    output_stream_name, field_name))

    ! check if we can handle it
    valid_input(b) = check_real_element_dim(domain % blocklist % allFields, &
      field_name, elem_name)

    if (valid_input(b)) then
      regions % number_of_variables = regions % number_of_variables + 1
    end if

    b = b + 1
  end do

  if (regions % number_of_variables < 1) then
    call mpas_dmpar_global_abort( &
      trim(CURRENT_CORE_NAME) // ' ERROR: there are no fields ' // &
      'in the regional stats output stream "' // trim(output_stream_name) // &
      '" that regional stats can be applied to.')
  end if

  ! create the memory
  allocate(regions % variables(regions % number_of_variables))

  ! create input variable name space
  do v = 1, regions % number_of_variables
    ! identifier
    write(var_identifier, '(I0)') v

    ! create input name space
    field_name = trim(storage_prefix) // '_' // trim(var_identifier) // &
      trim(INPUT_NAME_SUFFIX)
    call add_new_string(domain % blocklist % allFields, amPool, amPool, &
      field_name, regions % variables(v) % input_name)

    ! create has_vertical space
    field_name = trim(storage_prefix) // '_' // trim(var_identifier) // &
      trim(HAS_VERTICAL_SUFFIX)
    call add_new_integer(domain % blocklist % allFields, amPool, amPool, &
      field_name, regions % variables(v) % &
      has_vertical)
  end do

  ! get the old field names, assign to input name, and remove from stream
  call mpas_stream_mgr_begin_iteration(domain % streamManager, &
    output_stream_name, err)
  b = 1
  v = 1
  do while (mpas_stream_mgr_get_next_field(domain % streamManager, &
    output_stream_name, field_name))

    ! check if we can handle it
    if (valid_input(b)) then
      regions % variables(v) % input_name = field_name

      ! remove the old one from the stream
      call mpas_stream_mgr_remove_field(domain % streamManager, &
        output_stream_name, regions % variables(v) % input_name)

      v = v + 1
    end if

    b = b + 1
  end do

  deallocate(valid_input)

  !
  ! modify the stream
  !

  ! set up the variables and create output memory in stream
  call mpas_stream_mgr_begin_iteration(domain % streamManager, &
    output_stream_name, err)

  ! fetch region names and name via regions
  call mpas_pool_get_array(maskPool, REGION_NAMES_NAME, names)

  last = regions % num_regions_per(regions % group_index)
  ! allocate count and weight names
  allocate(regions % count_zerod_names(last))
  allocate(regions % weight_zerod_names(last))
  allocate(regions % count_oned_names(last))
  allocate(regions % weight_oned_names(last))
  do b = 1, last
    m = regions % groups(b, regions % group_index)

    ! counts
    regions % count_zerod_names(b) = output_naming(COUNT_TOKEN, &
      ZEROD_TOKEN, fix_region_name(names(m)), instance)

    call add_new_integer(domain % blocklist % allFields, amPool, amPool, &
      regions % count_zerod_names(b))

    ! add the field to the output stream
    call mpas_stream_mgr_add_field(domain % streamManager, &
      output_stream_name, regions % count_zerod_names(b), &
      ierr=err)

    ! create the counts and add to pool
    if (active_vertical) then
      ! count name
      regions % count_oned_names(b) = output_naming(COUNT_TOKEN, &
        ONED_TOKEN, fix_region_name(names(m)), instance)

      call add_new_integer_1d(domain % blocklist % allFields, &
        domain % blocklist % allFields, amPool, &
        regions % vertical_mask, regions % count_oned_names(b), &
        elem_name)

      ! add the field to the output stream
      call mpas_stream_mgr_add_field(domain % streamManager, &
        output_stream_name, regions % count_oned_names(b), &
        ierr=err)
    end if

    ! weights
    if (regions % function_oned == MUL_FUNC) then
      ! weight name
      regions % weight_zerod_names(b) = &
        output_naming(WEIGHT_TOKEN, &
        ZEROD_TOKEN, fix_region_name(names(m)), instance)

      call add_new_real(domain % blocklist % allFields, amPool, amPool, &
        regions % weight_zerod_names(b))

      ! add the field to the output stream
      call mpas_stream_mgr_add_field(domain % streamManager, &
        output_stream_name, regions % weight_zerod_names(b), &
        ierr=err)
    end if

    if (regions % function_twod == MUL_FUNC) then
      ! weight name
      regions % weight_oned_names(b) = &
        output_naming(WEIGHT_TOKEN, &
        ONED_TOKEN, fix_region_name(names(m)), instance)

      call add_new_real_1d(domain % blocklist % allFields, &
        domain % blocklist % allFields, amPool, &
        regions % vertical_mask, regions % weight_oned_names(b), &
        elem_name)

      ! add the field to the output stream
      call mpas_stream_mgr_add_field(domain % streamManager, &
        output_stream_name, regions % weight_oned_names(b), &
        ierr=err)
    end if
  end do

  do v = 1, regions % number_of_variables
    ! allocate output names
    allocate(regions % variables(v) % output_names(last))

    ! get the info of the field
    call mpas_pool_get_field_info(domain % blocklist % allFields, &
      regions % variables(v) % input_name, info)

    ! allocate output fields if the mask is active for this group
    do b = 1, last
      m = regions % groups(b, regions % group_index)

      ! region name
      regions % variables(v) % output_names(b) = output_naming(op_name, &
        regions % variables(v) % input_name, &
        fix_region_name(names(m)), instance)

      ! create the field and add to pool,
      ! also, set has_vertical if necessary
      if (b == 1) then
        if (active_vertical) then
          call add_new_field(info, domain % blocklist % allFields, &
            domain % blocklist % allFields, amPool, &
            regions % variables(v) % input_name, &
            regions % variables(v) % output_names(b), elem_name, &
            regions % variables(v) % has_vertical, regions % vertical_dim)
        else
          call add_new_field(info, domain % blocklist % allFields, &
            domain % blocklist % allFields, amPool, &
            regions % variables(v) % input_name, &
            regions % variables(v) % output_names(b), elem_name)

          regions % variables(v) % has_vertical = 0
        end if
      else
        call add_new_field(info, domain % blocklist % allFields, &
          domain % blocklist % allFields, amPool, &
          regions % variables(v) % input_name, &
          regions % variables(v) % output_names(b), elem_name)
      end if

      ! add the field to the output stream
      call mpas_stream_mgr_add_field(domain % streamManager, &
        output_stream_name, regions % variables(v) % output_names(b), &
        ierr=err)
    end do
  end do ! number_of_variables

  ! clean up the instance memory
  do v = 1, regions % number_of_variables
    deallocate(regions % variables(v) % output_names)
  end do
  deallocate(regions % variables)
  deallocate(regions % count_zerod_names)
  deallocate(regions % weight_zerod_names)
  deallocate(regions % count_oned_names)
  deallocate(regions % weight_oned_names)
end subroutine start_state



!***********************************************************************
! routine add_new_field
!
!> \brief Function to create a new field from an existing field
!> \author  Jon Woodring
!> \date    December 17, 2015
!> \details
!>  This routine conducts all initializations required for
!>  duplicating a field and adding it to a pool.
!-----------------------------------------------------------------------
subroutine add_new_field(info, all_fields, inpool, outpool, &
    inname, outname, elem_name, &
    has_vertical, vertical_dim)!{{{
  ! input variables
  type (mpas_pool_field_info_type), intent(in) :: info
  type (mpas_pool_type), pointer, intent(inout) :: inpool, outpool, all_fields
  character (len=StrKIND), intent(in) :: inname, outname, elem_name
  integer, intent(out), optional :: has_vertical
  character (len=StrKIND), intent(in), optional :: vertical_dim

  ! input/output variables

  ! output variables

  ! local variables

  ! duplicate field and add new field to inpool, outpool
  if ((info % nDims == 0) .or. (info % fieldType /= MPAS_POOL_REAL)) then
    call mpas_dmpar_global_abort( &
      trim(CURRENT_CORE_NAME) // ' ERROR: the impossible happened. ' // &
      'Tried to create an output field for an input field that ' // &
      'regional stats AM cannot support.')
  end if

  if (info % nDims == 1) then
    if (present(has_vertical)) then
      call copy_field_1r(all_fields, inpool, outpool, &
        inname, outname, elem_name, &
        has_vertical, vertical_dim)
    else
      call copy_field_1r(all_fields, inpool, outpool, &
        inname, outname, elem_name)
    end if
  else if (info % nDims == 2) then
    if (present(has_vertical)) then
      call copy_field_2r(all_fields, inpool, outpool, &
        inname, outname, elem_name, &
        has_vertical, vertical_dim)
    else
      call copy_field_2r(all_fields, inpool, outpool, &
        inname, outname, elem_name)
    end if
  else if (info % nDims == 3) then
    if (present(has_vertical)) then
      call copy_field_3r(all_fields, inpool, outpool, &
        inname, outname, elem_name, &
        has_vertical, vertical_dim)
    else
      call copy_field_3r(all_fields, inpool, outpool, &
        inname, outname, elem_name)
    end if
  else if (info % nDims == 4) then
    if (present(has_vertical)) then
      call copy_field_4r(all_fields, inpool, outpool, &
        inname, outname, elem_name, &
        has_vertical, vertical_dim)
    else
      call copy_field_4r(all_fields, inpool, outpool, &
        inname, outname, elem_name)
    end if
  else
    if (present(has_vertical)) then
      call copy_field_5r(all_fields, inpool, outpool, &
        inname, outname, elem_name, &
        has_vertical, vertical_dim)
    else
      call copy_field_5r(all_fields, inpool, outpool, &
        inname, outname, elem_name)
    end if
  end if

end subroutine add_new_field!}}}



!***********************************************************************
! routine typed_operate
!
!> \brief Do the operation, but switch on run-time type
!> \author  Jon Woodring
!> \date    December 17, 2015
!> \details
!>  Since we don't know the type of the array, we need to do some
!>  run-time type switching based on the type of the array.
!-----------------------------------------------------------------------
subroutine typed_operate(dminfo, blocklist, regions, variable) !{{{

  ! input variables
  type (dm_info), pointer, intent(in) :: dminfo
  type (block_type), pointer, intent(in) :: blocklist
  type (regional_type), intent(in) :: regions
  type (regional_variable_type), intent(in) :: variable

  ! input/output variables

  ! output variables

  ! local variables
  type (mpas_pool_field_info_type) :: info
  integer, pointer :: levels

  ! get the info
  call mpas_pool_get_field_info(blocklist % allFields, &
    variable % input_name, info)

  if (variable % has_vertical /= 0) then
    call mpas_pool_get_dimension(blocklist % dimensions, &
      regions % vertical_dim, levels)
  end if

  ! switch based on the type, dimensionality, operation, and vertical dim
  if (info % nDims == 1) then
    if (regions % operation == AVG_OP) then
      call operate1r_avg(dminfo, blocklist, regions, variable)
    else if (regions % operation == MIN_OP) then
      call operate1r_min(dminfo, blocklist, regions, variable)
    else if (regions % operation == MAX_OP) then
      call operate1r_max(dminfo, blocklist, regions, variable)
    else if (regions % operation == SUM_OP) then
      call operate1r_sum(dminfo, blocklist, regions, variable)
    else if (regions % operation == SOS_OP) then
      call operate1r_sos(dminfo, blocklist, regions, variable)
    else
      call mpas_dmpar_global_abort(trim(CURRENT_CORE_NAME) // ' ERROR: ' // &
        'the impossible happened - tried to operate with an unknown ' // &
        'operator in the regional stats AM')
    end if
  else if (info % nDims == 2) then
    if (regions % operation == AVG_OP) then
      if (variable % has_vertical /= 0) then
        call operatevert2r_avg(dminfo, blocklist, regions, variable, levels)
      else
        call operate2r_avg(dminfo, blocklist, regions, variable)
      end if
    else if (regions % operation == MIN_OP) then
      if (variable % has_vertical /= 0) then
        call operatevert2r_min(dminfo, blocklist, regions, variable, levels)
      else
        call operate2r_min(dminfo, blocklist, regions, variable)
      end if
    else if (regions % operation == MAX_OP) then
      if (variable % has_vertical /= 0) then
        call operatevert2r_max(dminfo, blocklist, regions, variable, levels)
      else
        call operate2r_max(dminfo, blocklist, regions, variable)
      end if
    else if (regions % operation == SUM_OP) then
      if (variable % has_vertical /= 0) then
        call operatevert2r_sum(dminfo, blocklist, regions, variable, levels)
      else
        call operate2r_sum(dminfo, blocklist, regions, variable)
      end if
    else if (regions % operation == SOS_OP) then
      if (variable % has_vertical /= 0) then
        call operatevert2r_sos(dminfo, blocklist, regions, variable, levels)
      else
        call operate2r_sos(dminfo, blocklist, regions, variable)
      end if
    else
      call mpas_dmpar_global_abort(trim(CURRENT_CORE_NAME) // ' ERROR: ' // &
        'the impossible happened - tried to operate with an unknown ' // &
        'operator in the regional stats AM')
    end if
  else if (info % nDims == 3) then
    if (regions % operation == AVG_OP) then
      if (variable % has_vertical /= 0) then
        call operatevert3r_avg(dminfo, blocklist, regions, variable, levels)
      else
        call operate3r_avg(dminfo, blocklist, regions, variable)
      end if
    else if (regions % operation == MIN_OP) then
      if (variable % has_vertical /= 0) then
        call operatevert3r_min(dminfo, blocklist, regions, variable, levels)
      else
        call operate3r_min(dminfo, blocklist, regions, variable)
      end if
    else if (regions % operation == MAX_OP) then
      if (variable % has_vertical /= 0) then
        call operatevert3r_max(dminfo, blocklist, regions, variable, levels)
      else
        call operate3r_max(dminfo, blocklist, regions, variable)
      end if
    else if (regions % operation == SUM_OP) then
      if (variable % has_vertical /= 0) then
        call operatevert3r_sum(dminfo, blocklist, regions, variable, levels)
      else
        call operate3r_sum(dminfo, blocklist, regions, variable)
      end if
    else if (regions % operation == SOS_OP) then
      if (variable % has_vertical /= 0) then
        call operatevert3r_sos(dminfo, blocklist, regions, variable, levels)
      else
        call operate3r_sos(dminfo, blocklist, regions, variable)
      end if
    else
      call mpas_dmpar_global_abort(trim(CURRENT_CORE_NAME) // ' ERROR: ' // &
        'the impossible happened - tried to operate with an unknown ' // &
        'operator in the regional stats AM')
    end if
  else if (info % nDims == 4) then
    if (regions % operation == AVG_OP) then
      if (variable % has_vertical /= 0) then
        call operatevert4r_avg(dminfo, blocklist, regions, variable, levels)
      else
        call operate4r_avg(dminfo, blocklist, regions, variable)
      end if
    else if (regions % operation == MIN_OP) then
      if (variable % has_vertical /= 0) then
        call operatevert4r_min(dminfo, blocklist, regions, variable, levels)
      else
        call operate4r_min(dminfo, blocklist, regions, variable)
      end if
    else if (regions % operation == MAX_OP) then
      if (variable % has_vertical /= 0) then
        call operatevert4r_max(dminfo, blocklist, regions, variable, levels)
      else
        call operate4r_max(dminfo, blocklist, regions, variable)
      end if
    else if (regions % operation == SUM_OP) then
      if (variable % has_vertical /= 0) then
        call operatevert4r_sum(dminfo, blocklist, regions, variable, levels)
      else
        call operate4r_sum(dminfo, blocklist, regions, variable)
      end if
    else if (regions % operation == SOS_OP) then
      if (variable % has_vertical /= 0) then
        call operatevert4r_sos(dminfo, blocklist, regions, variable, levels)
      else
        call operate4r_sos(dminfo, blocklist, regions, variable)
      end if
    else
      call mpas_dmpar_global_abort(trim(CURRENT_CORE_NAME) // ' ERROR: ' // &
        'the impossible happened - tried to operate with an unknown ' // &
        'operator in the regional stats AM')
    end if
  else if (info % nDims == 5) then
    if (regions % operation == AVG_OP) then
      if (variable % has_vertical /= 0) then
        call operatevert5r_avg(dminfo, blocklist, regions, variable, levels)
      else
        call operate5r_avg(dminfo, blocklist, regions, variable)
      end if
    else if (regions % operation == MIN_OP) then
      if (variable % has_vertical /= 0) then
        call operatevert5r_min(dminfo, blocklist, regions, variable, levels)
      else
        call operate5r_min(dminfo, blocklist, regions, variable)
      end if
    else if (regions % operation == MAX_OP) then
      if (variable % has_vertical /= 0) then
        call operatevert5r_max(dminfo, blocklist, regions, variable, levels)
      else
        call operate5r_max(dminfo, blocklist, regions, variable)
      end if
    else if (regions % operation == SUM_OP) then
      if (variable % has_vertical /= 0) then
        call operatevert5r_sum(dminfo, blocklist, regions, variable, levels)
      else
        call operate5r_sum(dminfo, blocklist, regions, variable)
      end if
    else if (regions % operation == SOS_OP) then
      if (variable % has_vertical /= 0) then
        call operatevert5r_sos(dminfo, blocklist, regions, variable, levels)
      else
        call operate5r_sos(dminfo, blocklist, regions, variable)
      end if
    else
      call mpas_dmpar_global_abort(trim(CURRENT_CORE_NAME) // ' ERROR: ' // &
        'the impossible happened - tried to operate with an unknown ' // &
        'operator in the regional stats AM')
    end if
  else
    call mpas_dmpar_global_abort(trim(CURRENT_CORE_NAME) // ' ERROR: ' // &
      'the impossible happened - tried to operate on a real field "' // &
      trim(variable % input_name) // '"that does not have 1-5 ' // &
      'dimensionality in the regional stats AM')
  end if

end subroutine typed_operate!}}}



!***********************************************************************
! routine copy_field_X
!
!> \brief Functions to create a new N-1 D field from an existing ND field
!> \author  Jon Woodring
!> \date    December 17, 2015
!> \details
!>  This routine conducts initializations required for
!>  duplicating a ND field and adding a N-1D to the AM's pool.
!>  This will only create an output in the first block of a blocklist.
!-----------------------------------------------------------------------
subroutine copy_field_1r(all_fields, inpool, outpool, &
    inname, outname, elem_name, &
    has_vertical, vertical_dim)!{{{
#include "regional_stats_inc/regional_field_1d_1.inc"
  type (field1DReal), pointer :: src
  type (field0DReal), pointer :: dst
#include "regional_stats_inc/regional_field_1d_2.inc"
end subroutine copy_field_1r!}}}


subroutine copy_field_2r(all_fields, inpool, outpool, &
    inname, outname, elem_name, &
    has_vertical, vertical_dim)!{{{
#include "regional_stats_inc/regional_field_nd_1.inc"
  type (field2DReal), pointer :: src
  type (field1DReal), pointer :: dst
  integer, dimension(2) :: src_dims
#include "regional_stats_inc/regional_field_nd_2.inc"
  allocate(dst % array(src_dims(1)))
#include "regional_stats_inc/regional_field_nd_3.inc"
end subroutine copy_field_2r!}}}


subroutine copy_field_3r(all_fields, inpool, outpool, &
    inname, outname, elem_name, &
    has_vertical, vertical_dim)!{{{
#include "regional_stats_inc/regional_field_nd_1.inc"
  type (field3DReal), pointer :: src
  type (field2DReal), pointer :: dst
  integer, dimension(3) :: src_dims
#include "regional_stats_inc/regional_field_nd_2.inc"
  allocate(dst % array(src_dims(1), src_dims(2)))
#include "regional_stats_inc/regional_field_nd_3.inc"
end subroutine copy_field_3r!}}}


subroutine copy_field_4r(all_fields, inpool, outpool, &
    inname, outname, elem_name, &
    has_vertical, vertical_dim)!{{{
#include "regional_stats_inc/regional_field_nd_1.inc"
  type (field4DReal), pointer :: src
  type (field3DReal), pointer :: dst
  integer, dimension(4) :: src_dims
#include "regional_stats_inc/regional_field_nd_2.inc"
  allocate(dst % array(src_dims(1), src_dims(2), src_dims(3)))
#include "regional_stats_inc/regional_field_nd_3.inc"
end subroutine copy_field_4r!}}}


subroutine copy_field_5r(all_fields, inpool, outpool, &
    inname, outname, elem_name, &
    has_vertical, vertical_dim)!{{{
#include "regional_stats_inc/regional_field_nd_1.inc"
  type (field5DReal), pointer :: src
  type (field4DReal), pointer :: dst
  integer, dimension(5) :: src_dims
#include "regional_stats_inc/regional_field_nd_2.inc"
  allocate(dst % array(src_dims(1), src_dims(2), src_dims(3), src_dims(4)))
#include "regional_stats_inc/regional_field_nd_3.inc"
end subroutine copy_field_5r!}}}



!***********************************************************************
! routine operateX_Y
!
!> \brief Series of subroutines to support operations on run-time types
!> \author  Jon Woodring
!> \date    December 17, 2015
!> \details
!>  These subroutines encapsulate the different operations that can occur
!>  based on the run-time types. (This would likely be
!>  instantiated generics/templates in other languages.)
!>
!>  We will do the reductions using slicing operations rather
!>  than for loops, and then copy to and from a flattened array,
!>  during distributed communication.
!>
!>  This will be slower than for loops using flattened array output
!>  indexing on the first copy, because we're going to have to copy twice
!>  on the distributed communication. (Because the framework reduceAll
!>  only works for rank-1 arrays.) The reason to do this is that
!>  to preserve the dimensionality semantics for legibility and
!>  for fewer bugs in the reduction. Basically, we flatten and unflatten,
!>  and could remove the flatten, and reduce in a for loop with indexing
!>  on a flattened array, and then reshape. While it could be possible to
!>  do manual for loops without slicing, and increasing the code,
!>  I doubt the performance gain for doing one less copy of a small
!>  array will be that significant.
!-----------------------------------------------------------------------

subroutine operate2r_avg (dminfo, start_block, regions, variable)
#include "regional_stats_inc/regional_op_start_1.inc"
  real (kind=RKIND), dimension(:, :), pointer :: in_array
  real (kind=RKIND), dimension(:), pointer :: out_array
#include "regional_stats_inc/regional_op_start_2.inc"
  out_array = 0
#include "regional_stats_inc/regional_op_start_3.inc"

#include "regional_stats_inc/regional_op_avg_1.inc"
  out_array = out_array + in_array(:, i) * mask(m, i)
#include "regional_stats_inc/regional_op_avg_2.inc"
  out_array = out_array + in_array(:, i) * (weights(i) * mask(m, i))
#include "regional_stats_inc/regional_op_avg_3.inc"
end subroutine operate2r_avg


subroutine operate5r_avg (dminfo, start_block, regions, variable)
#include "regional_stats_inc/regional_op_start_1.inc"
  real (kind=RKIND), dimension(:, :, :, :, :), pointer :: in_array
  real (kind=RKIND), dimension(:, :, :, :), pointer :: out_array
#include "regional_stats_inc/regional_op_start_2.inc"
  out_array = 0
#include "regional_stats_inc/regional_op_start_3.inc"

#include "regional_stats_inc/regional_op_avg_1.inc"
  out_array = out_array + in_array(:, :, :, :, i) * mask(m, i)
#include "regional_stats_inc/regional_op_avg_2.inc"
  out_array = out_array + in_array(:, :, :, :, i) * (weights(i) * mask(m, i))
#include "regional_stats_inc/regional_op_avg_3.inc"
end subroutine operate5r_avg


subroutine operate4r_avg (dminfo, start_block, regions, variable)
#include "regional_stats_inc/regional_op_start_1.inc"
  real (kind=RKIND), dimension(:, :, :, :), pointer :: in_array
  real (kind=RKIND), dimension(:, :, :), pointer :: out_array
#include "regional_stats_inc/regional_op_start_2.inc"
  out_array = 0
#include "regional_stats_inc/regional_op_start_3.inc"

#include "regional_stats_inc/regional_op_avg_1.inc"
  out_array = out_array + in_array(:, :, :, i) * mask(m, i)
#include "regional_stats_inc/regional_op_avg_2.inc"
  out_array = out_array + in_array(:, :, :, i) * (weights(i) * mask(m, i))
#include "regional_stats_inc/regional_op_avg_3.inc"
end subroutine operate4r_avg


subroutine operate3r_avg (dminfo, start_block, regions, variable)
#include "regional_stats_inc/regional_op_start_1.inc"
  real (kind=RKIND), dimension(:, :, :), pointer :: in_array
  real (kind=RKIND), dimension(:, :), pointer :: out_array
#include "regional_stats_inc/regional_op_start_2.inc"
  out_array = 0
#include "regional_stats_inc/regional_op_start_3.inc"

#include "regional_stats_inc/regional_op_avg_1.inc"
  out_array = out_array + in_array(:, :, i) * mask(m, i)
#include "regional_stats_inc/regional_op_avg_2.inc"
  out_array = out_array + in_array(:, :, i) * (weights(i) * mask(m, i))
#include "regional_stats_inc/regional_op_avg_3.inc"
end subroutine operate3r_avg


subroutine operate1r_avg (dminfo, start_block, regions, variable)
#include "regional_stats_inc/regional_op_start_1.inc"
  real (kind=RKIND), dimension(:), pointer :: in_array
  real (kind=RKIND), pointer :: out_array
#include "regional_stats_inc/regional_op_start_2.inc"
  out_array = 0
#include "regional_stats_inc/regional_op_start_3.inc"

#include "regional_stats_inc/regional_op_avg1d_1.inc"
  out_array = out_array + in_array(i) * mask(m, i)
#include "regional_stats_inc/regional_op_avg1d_2.inc"
  out_array = out_array + in_array(i) * (weights(i) * mask(m, i))
#include "regional_stats_inc/regional_op_avg1d_3.inc"
end subroutine operate1r_avg


subroutine operate2r_min (dminfo, start_block, regions, variable)
#include "regional_stats_inc/regional_op_start_1.inc"
  real (kind=RKIND), dimension(:, :), pointer :: in_array
  real (kind=RKIND), dimension(:), pointer :: out_array
#include "regional_stats_inc/regional_op_start_2.inc"
  out_array = DEFAULT_MPAS_MAX_VALUE
#include "regional_stats_inc/regional_op_start_3.inc"

#include "regional_stats_inc/regional_op_min_1.inc"
  out_array = min(out_array, &
    in_array(:, i) * mask(m, i) + &
    DEFAULT_MPAS_MAX_VALUE * (1 - mask(m, i)))
#include "regional_stats_inc/regional_op_min_2.inc"
end subroutine operate2r_min


subroutine operate5r_min (dminfo, start_block, regions, variable)
#include "regional_stats_inc/regional_op_start_1.inc"
  real (kind=RKIND), dimension(:, :, :, :, :), pointer :: in_array
  real (kind=RKIND), dimension(:, :, :, :), pointer :: out_array
#include "regional_stats_inc/regional_op_start_2.inc"
  out_array = DEFAULT_MPAS_MAX_VALUE
#include "regional_stats_inc/regional_op_start_3.inc"

#include "regional_stats_inc/regional_op_min_1.inc"
  out_array = min(out_array, &
    in_array(:, :, :, :, i) * mask(m, i) + &
    DEFAULT_MPAS_MAX_VALUE * (1 - mask(m, i)))
#include "regional_stats_inc/regional_op_min_2.inc"
end subroutine operate5r_min


subroutine operate4r_min (dminfo, start_block, regions, variable)
#include "regional_stats_inc/regional_op_start_1.inc"
  real (kind=RKIND), dimension(:, :, :, :), pointer :: in_array
  real (kind=RKIND), dimension(:, :, :), pointer :: out_array
#include "regional_stats_inc/regional_op_start_2.inc"
  out_array = DEFAULT_MPAS_MAX_VALUE
#include "regional_stats_inc/regional_op_start_3.inc"

#include "regional_stats_inc/regional_op_min_1.inc"
  out_array = min(out_array, &
    in_array(:, :, :, i) * mask(m, i) + &
    DEFAULT_MPAS_MAX_VALUE * (1 - mask(m, i)))
#include "regional_stats_inc/regional_op_min_2.inc"
end subroutine operate4r_min


subroutine operate3r_min (dminfo, start_block, regions, variable)
#include "regional_stats_inc/regional_op_start_1.inc"
  real (kind=RKIND), dimension(:, :, :), pointer :: in_array
  real (kind=RKIND), dimension(:, :), pointer :: out_array
#include "regional_stats_inc/regional_op_start_2.inc"
  out_array = DEFAULT_MPAS_MAX_VALUE
#include "regional_stats_inc/regional_op_start_3.inc"

#include "regional_stats_inc/regional_op_min_1.inc"
  out_array = min(out_array, &
    in_array(:, :, i) * mask(m, i) + &
    DEFAULT_MPAS_MAX_VALUE * (1 - mask(m, i)))
#include "regional_stats_inc/regional_op_min_2.inc"
end subroutine operate3r_min


subroutine operate1r_min (dminfo, start_block, regions, variable)
#include "regional_stats_inc/regional_op_start_1.inc"
  real (kind=RKIND), dimension(:), pointer :: in_array
  real (kind=RKIND), pointer :: out_array
#include "regional_stats_inc/regional_op_start_2.inc"
  out_array = DEFAULT_MPAS_MAX_VALUE
#include "regional_stats_inc/regional_op_start_3.inc"

#include "regional_stats_inc/regional_op_min1d_1.inc"
  out_array = min(out_array, &
    in_array(i) * mask(m, i) + &
    DEFAULT_MPAS_MAX_VALUE * (1 - mask(m, i)))
#include "regional_stats_inc/regional_op_min1d_2.inc"
end subroutine operate1r_min


subroutine operate2r_max (dminfo, start_block, regions, variable)
#include "regional_stats_inc/regional_op_start_1.inc"
  real (kind=RKIND), dimension(:, :), pointer :: in_array
  real (kind=RKIND), dimension(:), pointer :: out_array
#include "regional_stats_inc/regional_op_start_2.inc"
  out_array = DEFAULT_MPAS_MIN_VALUE
#include "regional_stats_inc/regional_op_start_3.inc"

#include "regional_stats_inc/regional_op_max_1.inc"
  out_array = max(out_array, &
    in_array(:, i) * mask(m, i) + &
    DEFAULT_MPAS_MIN_VALUE * (1 - mask(m, i)))
#include "regional_stats_inc/regional_op_max_2.inc"
end subroutine operate2r_max


subroutine operate5r_max (dminfo, start_block, regions, variable)
#include "regional_stats_inc/regional_op_start_1.inc"
  real (kind=RKIND), dimension(:, :, :, :, :), pointer :: in_array
  real (kind=RKIND), dimension(:, :, :, :), pointer :: out_array
#include "regional_stats_inc/regional_op_start_2.inc"
  out_array = DEFAULT_MPAS_MIN_VALUE
#include "regional_stats_inc/regional_op_start_3.inc"

#include "regional_stats_inc/regional_op_max_1.inc"
  out_array = max(out_array, &
    in_array(:, :, :, :, i) * mask(m, i) + &
    DEFAULT_MPAS_MIN_VALUE * (1 - mask(m, i)))
#include "regional_stats_inc/regional_op_max_2.inc"
end subroutine operate5r_max


subroutine operate4r_max (dminfo, start_block, regions, variable)
#include "regional_stats_inc/regional_op_start_1.inc"
  real (kind=RKIND), dimension(:, :, :, :), pointer :: in_array
  real (kind=RKIND), dimension(:, :, :), pointer :: out_array
#include "regional_stats_inc/regional_op_start_2.inc"
  out_array = DEFAULT_MPAS_MIN_VALUE
#include "regional_stats_inc/regional_op_start_3.inc"

#include "regional_stats_inc/regional_op_max_1.inc"
  out_array = max(out_array, &
    in_array(:, :, :, i) * mask(m, i) + &
    DEFAULT_MPAS_MIN_VALUE * (1 - mask(m, i)))
#include "regional_stats_inc/regional_op_max_2.inc"
end subroutine operate4r_max


subroutine operate3r_max (dminfo, start_block, regions, variable)
#include "regional_stats_inc/regional_op_start_1.inc"
  real (kind=RKIND), dimension(:, :, :), pointer :: in_array
  real (kind=RKIND), dimension(:, :), pointer :: out_array
#include "regional_stats_inc/regional_op_start_2.inc"
  out_array = DEFAULT_MPAS_MIN_VALUE
#include "regional_stats_inc/regional_op_start_3.inc"

#include "regional_stats_inc/regional_op_max_1.inc"
  out_array = max(out_array, &
    in_array(:, :, i) * mask(m, i) + &
    DEFAULT_MPAS_MIN_VALUE * (1 - mask(m, i)))
#include "regional_stats_inc/regional_op_max_2.inc"
end subroutine operate3r_max


subroutine operate1r_max (dminfo, start_block, regions, variable)
#include "regional_stats_inc/regional_op_start_1.inc"
  real (kind=RKIND), dimension(:), pointer :: in_array
  real (kind=RKIND), pointer :: out_array
#include "regional_stats_inc/regional_op_start_2.inc"
  out_array = DEFAULT_MPAS_MIN_VALUE
#include "regional_stats_inc/regional_op_start_3.inc"

#include "regional_stats_inc/regional_op_max1d_1.inc"
  out_array = max(out_array, &
    in_array(i) * mask(m, i) + &
    DEFAULT_MPAS_MIN_VALUE * (1 - mask(m, i)))
#include "regional_stats_inc/regional_op_max1d_2.inc"
end subroutine operate1r_max


subroutine operate2r_sum (dminfo, start_block, regions, variable)
#include "regional_stats_inc/regional_op_start_1.inc"
  real (kind=RKIND), dimension(:, :), pointer :: in_array
  real (kind=RKIND), dimension(:), pointer :: out_array
#include "regional_stats_inc/regional_op_start_2.inc"
  out_array = 0
#include "regional_stats_inc/regional_op_start_3.inc"

#include "regional_stats_inc/regional_op_sumsos_1.inc"
  out_array = out_array + in_array(:, i) * mask(m, i)
#include "regional_stats_inc/regional_op_sumsos_2.inc"
  out_array = out_array + in_array(:, i) * (weights(i) * mask(m, i))
#include "regional_stats_inc/regional_op_sumsos_3.inc"
end subroutine operate2r_sum


subroutine operate5r_sum (dminfo, start_block, regions, variable)
#include "regional_stats_inc/regional_op_start_1.inc"
  real (kind=RKIND), dimension(:, :, :, :, :), pointer :: in_array
  real (kind=RKIND), dimension(:, :, :, :), pointer :: out_array
#include "regional_stats_inc/regional_op_start_2.inc"
  out_array = 0
#include "regional_stats_inc/regional_op_start_3.inc"

#include "regional_stats_inc/regional_op_sumsos_1.inc"
  out_array = out_array + in_array(:, :, :, :, i) * mask(m, i)
#include "regional_stats_inc/regional_op_sumsos_2.inc"
  out_array = out_array + in_array(:, :, :, :, i) * (weights(i) * mask(m, i))
#include "regional_stats_inc/regional_op_sumsos_3.inc"
end subroutine operate5r_sum


subroutine operate4r_sum (dminfo, start_block, regions, variable)
#include "regional_stats_inc/regional_op_start_1.inc"
  real (kind=RKIND), dimension(:, :, :, :), pointer :: in_array
  real (kind=RKIND), dimension(:, :, :), pointer :: out_array
#include "regional_stats_inc/regional_op_start_2.inc"
  out_array = 0
#include "regional_stats_inc/regional_op_start_3.inc"

#include "regional_stats_inc/regional_op_sumsos_1.inc"
  out_array = out_array + in_array(:, :, :, i) * mask(m, i)
#include "regional_stats_inc/regional_op_sumsos_2.inc"
  out_array = out_array + in_array(:, :, :, i) * (weights(i) * mask(m, i))
#include "regional_stats_inc/regional_op_sumsos_3.inc"
end subroutine operate4r_sum


subroutine operate3r_sum (dminfo, start_block, regions, variable)
#include "regional_stats_inc/regional_op_start_1.inc"
  real (kind=RKIND), dimension(:, :, :), pointer :: in_array
  real (kind=RKIND), dimension(:, :), pointer :: out_array
#include "regional_stats_inc/regional_op_start_2.inc"
  out_array = 0
#include "regional_stats_inc/regional_op_start_3.inc"

#include "regional_stats_inc/regional_op_sumsos_1.inc"
  out_array = out_array + in_array(:, :, i) * mask(m, i)
#include "regional_stats_inc/regional_op_sumsos_2.inc"
  out_array = out_array + in_array(:, :, i) * (weights(i) * mask(m, i))
#include "regional_stats_inc/regional_op_sumsos_3.inc"
end subroutine operate3r_sum


subroutine operate1r_sum (dminfo, start_block, regions, variable)
#include "regional_stats_inc/regional_op_start_1.inc"
  real (kind=RKIND), dimension(:), pointer :: in_array
  real (kind=RKIND), pointer :: out_array
#include "regional_stats_inc/regional_op_start_2.inc"
  out_array = 0
#include "regional_stats_inc/regional_op_start_3.inc"

#include "regional_stats_inc/regional_op_sumsos1d_1.inc"
  out_array = out_array + in_array(i) * mask(m, i)
#include "regional_stats_inc/regional_op_sumsos1d_2.inc"
  out_array = out_array + in_array(i) * (weights(i) * mask(m, i))
#include "regional_stats_inc/regional_op_sumsos1d_3.inc"
end subroutine operate1r_sum


subroutine operate2r_sos (dminfo, start_block, regions, variable)
#include "regional_stats_inc/regional_op_start_1.inc"
  real (kind=RKIND), dimension(:, :), pointer :: in_array
  real (kind=RKIND), dimension(:), pointer :: out_array
#include "regional_stats_inc/regional_op_start_2.inc"
  out_array = 0
#include "regional_stats_inc/regional_op_start_3.inc"

#include "regional_stats_inc/regional_op_sumsos_1.inc"
  out_array = out_array + in_array(:, i) * in_array(:, i) * mask(m, i)
#include "regional_stats_inc/regional_op_sumsos_2.inc"
  out_array = out_array + in_array(:, i) * &
    in_array(:, i) * (weights(i) * mask(m, i))
#include "regional_stats_inc/regional_op_sumsos_3.inc"
end subroutine operate2r_sos


subroutine operate5r_sos (dminfo, start_block, regions, variable)
#include "regional_stats_inc/regional_op_start_1.inc"
  real (kind=RKIND), dimension(:, :, :, :, :), pointer :: in_array
  real (kind=RKIND), dimension(:, :, :, :), pointer :: out_array
#include "regional_stats_inc/regional_op_start_2.inc"
  out_array = 0
#include "regional_stats_inc/regional_op_start_3.inc"

#include "regional_stats_inc/regional_op_sumsos_1.inc"
  out_array = out_array + in_array(:, :, :, :, i) * &
    in_array(:, :, :, :, i) * mask(m, i)
#include "regional_stats_inc/regional_op_sumsos_2.inc"
  out_array = out_array + in_array(:, :, :, :, i) * &
    in_array(:, :, :, :, i) * (weights(i) * mask(m, i))
#include "regional_stats_inc/regional_op_sumsos_3.inc"
end subroutine operate5r_sos


subroutine operate4r_sos (dminfo, start_block, regions, variable)
#include "regional_stats_inc/regional_op_start_1.inc"
  real (kind=RKIND), dimension(:, :, :, :), pointer :: in_array
  real (kind=RKIND), dimension(:, :, :), pointer :: out_array
#include "regional_stats_inc/regional_op_start_2.inc"
  out_array = 0
#include "regional_stats_inc/regional_op_start_3.inc"

#include "regional_stats_inc/regional_op_sumsos_1.inc"
  out_array = out_array + in_array(:, :, :, i) * &
    in_array(:, :, :, i) * mask(m, i)
#include "regional_stats_inc/regional_op_sumsos_2.inc"
  out_array = out_array + in_array(:, :, :, i) * &
    in_array(:, :, :, i) * (weights(i) * mask(m, i))
#include "regional_stats_inc/regional_op_sumsos_3.inc"
end subroutine operate4r_sos


subroutine operate3r_sos (dminfo, start_block, regions, variable)
#include "regional_stats_inc/regional_op_start_1.inc"
  real (kind=RKIND), dimension(:, :, :), pointer :: in_array
  real (kind=RKIND), dimension(:, :), pointer :: out_array
#include "regional_stats_inc/regional_op_start_2.inc"
  out_array = 0
#include "regional_stats_inc/regional_op_start_3.inc"

#include "regional_stats_inc/regional_op_sumsos_1.inc"
  out_array = out_array + in_array(:, :, i) * &
    in_array(:, :, i) * mask(m, i)
#include "regional_stats_inc/regional_op_sumsos_2.inc"
  out_array = out_array + in_array(:, :, i) * &
    in_array(:, :, i) * (weights(i) * mask(m, i))
#include "regional_stats_inc/regional_op_sumsos_3.inc"
end subroutine operate3r_sos


subroutine operate1r_sos (dminfo, start_block, regions, variable)
#include "regional_stats_inc/regional_op_start_1.inc"
  real (kind=RKIND), dimension(:), pointer :: in_array
  real (kind=RKIND), pointer :: out_array
#include "regional_stats_inc/regional_op_start_2.inc"
  out_array = 0
#include "regional_stats_inc/regional_op_start_3.inc"

#include "regional_stats_inc/regional_op_sumsos1d_1.inc"
  out_array = out_array + in_array(i) * &
    in_array(i) * mask(m, i)
#include "regional_stats_inc/regional_op_sumsos1d_2.inc"
  out_array = out_array + in_array(i) * &
    in_array(i) * (weights(i) * mask(m, i))
#include "regional_stats_inc/regional_op_sumsos1d_3.inc"
end subroutine operate1r_sos


subroutine operatevert2r_avg (dminfo, start_block, regions, variable, levels)
#include "regional_stats_inc/regional_opvert_start_1.inc"
  real (kind=RKIND), dimension(:, :), pointer :: in_array
  real (kind=RKIND), dimension(:), pointer :: out_array
#include "regional_stats_inc/regional_opvert_start_2.inc"
  out_array = 0
#include "regional_stats_inc/regional_opvert_start_3.inc"

#include "regional_stats_inc/regional_opvert_avg_1.inc"
  out_array(v) = out_array(v) + in_array(v, i) * &
    (mask(m, i) * vertical_mask(v, i))
#include "regional_stats_inc/regional_opvert_avg_2.inc"
  out_array(v) = out_array(v) + in_array(v, i) * &
    (weights(v, i) * (mask(m, i) * vertical_mask(v, i)))
#include "regional_stats_inc/regional_opvert_avg_3.inc"
  do v = 1, levels
    if (count_array(v) > 0) then
      if (regions % function_twod == ID_FUNC) then
        out_array(v) = out_array(v) / count_array(v)
      else
        out_array(v) = out_array(v) / weight_total(v)
      end if
    end if
  end do
#include "regional_stats_inc/regional_opvert_avg_4.inc"
end subroutine operatevert2r_avg


subroutine operatevert5r_avg (dminfo, start_block, regions, variable, levels)
#include "regional_stats_inc/regional_opvert_start_1.inc"
  real (kind=RKIND), dimension(:, :, :, :, :), pointer :: in_array
  real (kind=RKIND), dimension(:, :, :, :), pointer :: out_array
#include "regional_stats_inc/regional_opvert_start_2.inc"
  out_array = 0
#include "regional_stats_inc/regional_opvert_start_3.inc"

#include "regional_stats_inc/regional_opvert_avg_1.inc"
  out_array(:, :, :, v) = out_array(:, :, :, v) + in_array(:, :, :, v, i) * &
    (mask(m, i) * vertical_mask(v, i))
#include "regional_stats_inc/regional_opvert_avg_2.inc"
  out_array(:, :, :, v) = out_array(:, :, :, v) + in_array(:, :, :, v, i) * &
    (weights(v, i) * (mask(m, i) * vertical_mask(v, i)))
#include "regional_stats_inc/regional_opvert_avg_3.inc"
  do v = 1, levels
    if (count_array(v) > 0) then
      if (regions % function_twod == ID_FUNC) then
        out_array(:, :, :, v) = out_array(:, :, :, v) / count_array(v)
      else
        out_array(:, :, :, v) = out_array(:, :, :, v) / weight_total(v)
      end if
    end if
  end do
#include "regional_stats_inc/regional_opvert_avg_4.inc"
end subroutine operatevert5r_avg


subroutine operatevert4r_avg (dminfo, start_block, regions, variable, levels)
#include "regional_stats_inc/regional_opvert_start_1.inc"
  real (kind=RKIND), dimension(:, :, :, :), pointer :: in_array
  real (kind=RKIND), dimension(:, :, :), pointer :: out_array
#include "regional_stats_inc/regional_opvert_start_2.inc"
  out_array = 0
#include "regional_stats_inc/regional_opvert_start_3.inc"

#include "regional_stats_inc/regional_opvert_avg_1.inc"
  out_array(:, :, v) = out_array(:, :, v) + in_array(:, :, v, i) * &
    (mask(m, i) * vertical_mask(v, i))
#include "regional_stats_inc/regional_opvert_avg_2.inc"
  out_array(:, :, v) = out_array(:, :, v) + in_array(:, :, v, i) * &
    (weights(v, i) * (mask(m, i) * vertical_mask(v, i)))
#include "regional_stats_inc/regional_opvert_avg_3.inc"
  do v = 1, levels
    if (count_array(v) > 0) then
      if (regions % function_twod == ID_FUNC) then
        out_array(:, :, v) = out_array(:, :, v) / count_array(v)
      else
        out_array(:, :, v) = out_array(:, :, v) / weight_total(v)
      end if
    end if
  end do
#include "regional_stats_inc/regional_opvert_avg_4.inc"
end subroutine operatevert4r_avg


subroutine operatevert3r_avg (dminfo, start_block, regions, variable, levels)
#include "regional_stats_inc/regional_opvert_start_1.inc"
  real (kind=RKIND), dimension(:, :, :), pointer :: in_array
  real (kind=RKIND), dimension(:, :), pointer :: out_array
#include "regional_stats_inc/regional_opvert_start_2.inc"
  out_array = 0
#include "regional_stats_inc/regional_opvert_start_3.inc"

#include "regional_stats_inc/regional_opvert_avg_1.inc"
  out_array(:, v) = out_array(:, v) + in_array(:, v, i) * &
    (mask(m, i) * vertical_mask(v, i))
#include "regional_stats_inc/regional_opvert_avg_2.inc"
  out_array(:, v) = out_array(:, v) + in_array(:, v, i) * &
    (weights(v, i) * (mask(m, i) * vertical_mask(v, i)))
#include "regional_stats_inc/regional_opvert_avg_3.inc"
  do v = 1, levels
    if (count_array(v) > 0) then
      if (regions % function_twod == ID_FUNC) then
        out_array(:, v) = out_array(:, v) / count_array(v)
      else
        out_array(:, v) = out_array(:, v) / weight_total(v)
      end if
    end if
  end do
#include "regional_stats_inc/regional_opvert_avg_4.inc"
end subroutine operatevert3r_avg


subroutine operatevert2r_min (dminfo, start_block, regions, variable, levels)
#include "regional_stats_inc/regional_opvert_start_1.inc"
  real (kind=RKIND), dimension(:, :), pointer :: in_array
  real (kind=RKIND), dimension(:), pointer :: out_array
#include "regional_stats_inc/regional_opvert_start_2.inc"
  out_array = DEFAULT_MPAS_MAX_VALUE
#include "regional_stats_inc/regional_opvert_start_3.inc"

#include "regional_stats_inc/regional_opvert_min_1.inc"
  out_array(v) = min(out_array(v), &
    in_array(v, i) * (mask(m, i) * vertical_mask(v, i)) + &
    DEFAULT_MPAS_MAX_VALUE * (1 - (mask(m, i) * vertical_mask(v, i))))
#include "regional_stats_inc/regional_opvert_min_2.inc"
end subroutine operatevert2r_min


subroutine operatevert5r_min (dminfo, start_block, regions, variable, levels)
#include "regional_stats_inc/regional_opvert_start_1.inc"
  real (kind=RKIND), dimension(:, :, :, :, :), pointer :: in_array
  real (kind=RKIND), dimension(:, :, :, :), pointer :: out_array
#include "regional_stats_inc/regional_opvert_start_2.inc"
  out_array = DEFAULT_MPAS_MAX_VALUE
#include "regional_stats_inc/regional_opvert_start_3.inc"

#include "regional_stats_inc/regional_opvert_min_1.inc"
  out_array(:, :, :, v) = min(out_array(:, :, :, v), &
    in_array(:, :, :, v, i) * (mask(m, i) * vertical_mask(v, i)) + &
    DEFAULT_MPAS_MAX_VALUE * (1 - (mask(m, i) * vertical_mask(v, i))))
#include "regional_stats_inc/regional_opvert_min_2.inc"
end subroutine operatevert5r_min


subroutine operatevert4r_min (dminfo, start_block, regions, variable, levels)
#include "regional_stats_inc/regional_opvert_start_1.inc"
  real (kind=RKIND), dimension(:, :, :, :), pointer :: in_array
  real (kind=RKIND), dimension(:, :, :), pointer :: out_array
#include "regional_stats_inc/regional_opvert_start_2.inc"
  out_array = DEFAULT_MPAS_MAX_VALUE
#include "regional_stats_inc/regional_opvert_start_3.inc"

#include "regional_stats_inc/regional_opvert_min_1.inc"
  out_array(:, :, v) = min(out_array(:, :, v), &
    in_array(:, :, v, i) * (mask(m, i) * vertical_mask(v, i)) + &
    DEFAULT_MPAS_MAX_VALUE * (1 - (mask(m, i) * vertical_mask(v, i))))
#include "regional_stats_inc/regional_opvert_min_2.inc"
end subroutine operatevert4r_min


subroutine operatevert3r_min (dminfo, start_block, regions, variable, levels)
#include "regional_stats_inc/regional_opvert_start_1.inc"
  real (kind=RKIND), dimension(:, :, :), pointer :: in_array
  real (kind=RKIND), dimension(:, :), pointer :: out_array
#include "regional_stats_inc/regional_opvert_start_2.inc"
  out_array = DEFAULT_MPAS_MAX_VALUE
#include "regional_stats_inc/regional_opvert_start_3.inc"

#include "regional_stats_inc/regional_opvert_min_1.inc"
  out_array(:, v) = min(out_array(:, v), &
    in_array(:, v, i) * (mask(m, i) * vertical_mask(v, i)) + &
    DEFAULT_MPAS_MAX_VALUE * (1 - (mask(m, i) * vertical_mask(v, i))))
#include "regional_stats_inc/regional_opvert_min_2.inc"
end subroutine operatevert3r_min


subroutine operatevert2r_max (dminfo, start_block, regions, variable, levels)
#include "regional_stats_inc/regional_opvert_start_1.inc"
  real (kind=RKIND), dimension(:, :), pointer :: in_array
  real (kind=RKIND), dimension(:), pointer :: out_array
#include "regional_stats_inc/regional_opvert_start_2.inc"
  out_array = DEFAULT_MPAS_MIN_VALUE
#include "regional_stats_inc/regional_opvert_start_3.inc"

#include "regional_stats_inc/regional_opvert_max_1.inc"
  out_array(v) = max(out_array(v), &
    in_array(v, i) * (mask(m, i) * vertical_mask(v, i)) + &
    DEFAULT_MPAS_MIN_VALUE * (1 - (mask(m, i) * vertical_mask(v, i))))
#include "regional_stats_inc/regional_opvert_max_2.inc"
end subroutine operatevert2r_max


subroutine operatevert5r_max (dminfo, start_block, regions, variable, levels)
#include "regional_stats_inc/regional_opvert_start_1.inc"
  real (kind=RKIND), dimension(:, :, :, :, :), pointer :: in_array
  real (kind=RKIND), dimension(:, :, :, :), pointer :: out_array
#include "regional_stats_inc/regional_opvert_start_2.inc"
  out_array = DEFAULT_MPAS_MIN_VALUE
#include "regional_stats_inc/regional_opvert_start_3.inc"

#include "regional_stats_inc/regional_opvert_max_1.inc"
  out_array(:, :, :, v) = max(out_array(:, :, :, v), &
    in_array(:, :, :, v, i) * (mask(m, i) * vertical_mask(v, i)) + &
    DEFAULT_MPAS_MIN_VALUE * (1 - (mask(m, i) * vertical_mask(v, i))))
#include "regional_stats_inc/regional_opvert_max_2.inc"
end subroutine operatevert5r_max


subroutine operatevert4r_max (dminfo, start_block, regions, variable, levels)
#include "regional_stats_inc/regional_opvert_start_1.inc"
  real (kind=RKIND), dimension(:, :, :, :), pointer :: in_array
  real (kind=RKIND), dimension(:, :, :), pointer :: out_array
#include "regional_stats_inc/regional_opvert_start_2.inc"
  out_array = DEFAULT_MPAS_MIN_VALUE
#include "regional_stats_inc/regional_opvert_start_3.inc"

#include "regional_stats_inc/regional_opvert_max_1.inc"
  out_array(:, :, v) = max(out_array(:, :, v), &
    in_array(:, :, v, i) * (mask(m, i) * vertical_mask(v, i)) + &
    DEFAULT_MPAS_MIN_VALUE * (1 - (mask(m, i) * vertical_mask(v, i))))
#include "regional_stats_inc/regional_opvert_max_2.inc"
end subroutine operatevert4r_max


subroutine operatevert3r_max (dminfo, start_block, regions, variable, levels)
#include "regional_stats_inc/regional_opvert_start_1.inc"
  real (kind=RKIND), dimension(:, :, :), pointer :: in_array
  real (kind=RKIND), dimension(:, :), pointer :: out_array
#include "regional_stats_inc/regional_opvert_start_2.inc"
  out_array = DEFAULT_MPAS_MIN_VALUE
#include "regional_stats_inc/regional_opvert_start_3.inc"

#include "regional_stats_inc/regional_opvert_max_1.inc"
  out_array(:, v) = max(out_array(:, v), &
    in_array(:, v, i) * (mask(m, i) * vertical_mask(v, i)) + &
    DEFAULT_MPAS_MIN_VALUE * (1 - (mask(m, i) * vertical_mask(v, i))))
#include "regional_stats_inc/regional_opvert_max_2.inc"
end subroutine operatevert3r_max


subroutine operatevert2r_sum (dminfo, start_block, regions, variable, levels)
#include "regional_stats_inc/regional_opvert_start_1.inc"
  real (kind=RKIND), dimension(:, :), pointer :: in_array
  real (kind=RKIND), dimension(:), pointer :: out_array
#include "regional_stats_inc/regional_opvert_start_2.inc"
  out_array = 0
#include "regional_stats_inc/regional_opvert_start_3.inc"

#include "regional_stats_inc/regional_opvert_sumsos_1.inc"
  out_array(v) = out_array(v) + in_array(v, i) * &
    (mask(m, i) * vertical_mask(v, i))
#include "regional_stats_inc/regional_opvert_sumsos_2.inc"
  out_array(v) = out_array(v) + in_array(v, i) * &
    (weights(v, i) * (mask(m, i) * vertical_mask(v, i)))
#include "regional_stats_inc/regional_opvert_sumsos_3.inc"
end subroutine operatevert2r_sum


subroutine operatevert5r_sum (dminfo, start_block, regions, variable, levels)
#include "regional_stats_inc/regional_opvert_start_1.inc"
  real (kind=RKIND), dimension(:, :, :, :, :), pointer :: in_array
  real (kind=RKIND), dimension(:, :, :, :), pointer :: out_array
#include "regional_stats_inc/regional_opvert_start_2.inc"
  out_array = 0
#include "regional_stats_inc/regional_opvert_start_3.inc"

#include "regional_stats_inc/regional_opvert_sumsos_1.inc"
  out_array(:, :, :, v) = out_array(:, :, :, v) + in_array(:, :, :, v, i) * &
    (mask(m, i) * vertical_mask(v, i))
#include "regional_stats_inc/regional_opvert_sumsos_2.inc"
  out_array(:, :, :, v) = out_array(:, :, :, v) + in_array(:, :, :, v, i) * &
    (weights(v, i) * (mask(m, i) * vertical_mask(v, i)))
#include "regional_stats_inc/regional_opvert_sumsos_3.inc"
end subroutine operatevert5r_sum


subroutine operatevert4r_sum (dminfo, start_block, regions, variable, levels)
#include "regional_stats_inc/regional_opvert_start_1.inc"
  real (kind=RKIND), dimension(:, :, :, :), pointer :: in_array
  real (kind=RKIND), dimension(:, :, :), pointer :: out_array
#include "regional_stats_inc/regional_opvert_start_2.inc"
  out_array = 0
#include "regional_stats_inc/regional_opvert_start_3.inc"

#include "regional_stats_inc/regional_opvert_sumsos_1.inc"
  out_array(:, :, v) = out_array(:, :, v) + in_array(:, :, v, i) * &
    (mask(m, i) * vertical_mask(v, i))
#include "regional_stats_inc/regional_opvert_sumsos_2.inc"
  out_array(:, :, v) = out_array(:, :, v) + in_array(:, :, v, i) * &
    (weights(v, i) * (mask(m, i) * vertical_mask(v, i)))
#include "regional_stats_inc/regional_opvert_sumsos_3.inc"
end subroutine operatevert4r_sum


subroutine operatevert3r_sum (dminfo, start_block, regions, variable, levels)
#include "regional_stats_inc/regional_opvert_start_1.inc"
  real (kind=RKIND), dimension(:, :, :), pointer :: in_array
  real (kind=RKIND), dimension(:, :), pointer :: out_array
#include "regional_stats_inc/regional_opvert_start_2.inc"
  out_array = 0
#include "regional_stats_inc/regional_opvert_start_3.inc"

#include "regional_stats_inc/regional_opvert_sumsos_1.inc"
  out_array(:, v) = out_array(:, v) + in_array(:, v, i) * &
    (mask(m, i) * vertical_mask(v, i))
#include "regional_stats_inc/regional_opvert_sumsos_2.inc"
  out_array(:, v) = out_array(:, v) + in_array(:, v, i) * &
    (weights(v, i) * (mask(m, i) * vertical_mask(v, i)))
#include "regional_stats_inc/regional_opvert_sumsos_3.inc"
end subroutine operatevert3r_sum


subroutine operatevert2r_sos (dminfo, start_block, regions, variable, levels)
#include "regional_stats_inc/regional_opvert_start_1.inc"
  real (kind=RKIND), dimension(:, :), pointer :: in_array
  real (kind=RKIND), dimension(:), pointer :: out_array
#include "regional_stats_inc/regional_opvert_start_2.inc"
  out_array = 0
#include "regional_stats_inc/regional_opvert_start_3.inc"

#include "regional_stats_inc/regional_opvert_sumsos_1.inc"
  out_array(v) = out_array(v) + in_array(v, i) * &
    in_array(v, i) * (mask(m, i) * vertical_mask(v, i))
#include "regional_stats_inc/regional_opvert_sumsos_2.inc"
  out_array(v) = out_array(v) + in_array(v, i) * &
    in_array(v, i) * &
    (weights(v, i) * (mask(m, i) * vertical_mask(v, i)))
#include "regional_stats_inc/regional_opvert_sumsos_3.inc"
end subroutine operatevert2r_sos


subroutine operatevert5r_sos (dminfo, start_block, regions, variable, levels)
#include "regional_stats_inc/regional_opvert_start_1.inc"
  real (kind=RKIND), dimension(:, :, :, :, :), pointer :: in_array
  real (kind=RKIND), dimension(:, :, :, :), pointer :: out_array
#include "regional_stats_inc/regional_opvert_start_2.inc"
  out_array = 0
#include "regional_stats_inc/regional_opvert_start_3.inc"

#include "regional_stats_inc/regional_opvert_sumsos_1.inc"
  out_array(:, :, :, v) = out_array(:, :, :, v) + in_array(:, :, :, v, i) * &
    in_array(:, :, :, v, i) * (mask(m, i) * vertical_mask(v, i))
#include "regional_stats_inc/regional_opvert_sumsos_2.inc"
  out_array(:, :, :, v) = out_array(:, :, :, v) + in_array(:, :, :, v, i) * &
    in_array(:, :, :, v, i) * &
    (weights(v, i) * (mask(m, i) * vertical_mask(v, i)))
#include "regional_stats_inc/regional_opvert_sumsos_3.inc"
end subroutine operatevert5r_sos


subroutine operatevert4r_sos (dminfo, start_block, regions, variable, levels)
#include "regional_stats_inc/regional_opvert_start_1.inc"
  real (kind=RKIND), dimension(:, :, :, :), pointer :: in_array
  real (kind=RKIND), dimension(:, :, :), pointer :: out_array
#include "regional_stats_inc/regional_opvert_start_2.inc"
  out_array = 0
#include "regional_stats_inc/regional_opvert_start_3.inc"

#include "regional_stats_inc/regional_opvert_sumsos_1.inc"
  out_array(:, :, v) = out_array(:, :, v) + in_array(:, :, v, i) * &
    in_array(:, :, v, i) * (mask(m, i) * vertical_mask(v, i))
#include "regional_stats_inc/regional_opvert_sumsos_2.inc"
  out_array(:, :, v) = out_array(:, :, v) + in_array(:, :, v, i) * &
    in_array(:, :, v, i) * &
    (weights(v, i) * (mask(m, i) * vertical_mask(v, i)))
#include "regional_stats_inc/regional_opvert_sumsos_3.inc"
end subroutine operatevert4r_sos


subroutine operatevert3r_sos (dminfo, start_block, regions, variable, levels)
#include "regional_stats_inc/regional_opvert_start_1.inc"
  real (kind=RKIND), dimension(:, :, :), pointer :: in_array
  real (kind=RKIND), dimension(:, :), pointer :: out_array
#include "regional_stats_inc/regional_opvert_start_2.inc"
  out_array = 0
#include "regional_stats_inc/regional_opvert_start_3.inc"

#include "regional_stats_inc/regional_opvert_sumsos_1.inc"
  out_array(:, v) = out_array(:, v) + in_array(:, v, i) * &
    in_array(:, v, i) * (mask(m, i) * vertical_mask(v, i))
#include "regional_stats_inc/regional_opvert_sumsos_2.inc"
  out_array(:, v) = out_array(:, v) + in_array(:, v, i) * &
    in_array(:, v, i) * &
    (weights(v, i) * (mask(m, i) * vertical_mask(v, i)))
#include "regional_stats_inc/regional_opvert_sumsos_3.inc"
end subroutine operatevert3r_sos


end module ocn_regional_stats
! vim: foldmethod=marker
