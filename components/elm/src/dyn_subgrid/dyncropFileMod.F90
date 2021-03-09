module dyncropFileMod

  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Handle reading of the dataset that specifies transient areas the crop landunit as
  ! well as the breakdown of each crop.
  !
  ! !USES:
#include "shr_assert.h"
  use shr_log_mod           , only : errMsg => shr_log_errMsg
  use shr_kind_mod          , only : r8 => shr_kind_r8
  use decompMod             , only : bounds_type, BOUNDS_LEVEL_PROC
  use dynFileMod            , only : dyn_file_type
  use dynVarTimeUninterpMod , only : dyn_var_time_uninterp_type
  use elm_varctl            , only : iulog, use_crop
  use elm_varcon            , only : grlnd, namec
  use abortutils            , only : endrun
  use spmdMod               , only : masterproc, mpicom
  use LandunitType          , only : lun_pp
  use ColumnType            , only : col_pp
  use VegetationType        , only : veg_pp
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  private
  save
  public :: dyncrop_init     ! initialize information read from landuse.timeseries dataset
  public :: dyncrop_interp   ! get crop data for the current time step, if needed
  !
  ! ! PRIVATE TYPES
  type(dyn_file_type), target      :: dyncrop_file ! information for the file containing transient crop data
  type(dyn_var_time_uninterp_type) :: wtcrop       ! weight of the crop landunit
  type(dyn_var_time_uninterp_type) :: wtcft        ! weight of each CFT relative to the crop landunit
  type(dyn_var_time_uninterp_type) :: fertcft      ! fertilizer of each CFT

  ! Names of variables on file
  character(len=*), parameter :: crop_varname = 'PCT_CROP'
  character(len=*), parameter :: cft_varname  = 'PCT_CFT'
  character(len=*), parameter :: fert_varname  = 'FERTNITRO_CFT'

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !---------------------------------------------------------------------------

contains
  !-----------------------------------------------------------------------
  subroutine dyncrop_init(bounds, dyncrop_filename)
    !
    ! !DESCRIPTION:
    ! Initialize dataset containing transient crop info (position it to the right time
    ! samples that bound the initial model date)
    !
    ! !USES:
    use elm_varpar     , only : cft_size
    use ncdio_pio      , only : check_dim
    use dynTimeInfoMod , only : YEAR_POSITION_START_OF_TIMESTEP
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in) :: bounds           ! proc-level bounds
    character(len=*)  , intent(in) :: dyncrop_filename ! name of file containing transient crop information
    !
    ! !LOCAL VARIABLES:
    integer :: num_points     ! number of spatial points
    integer :: wtcft_shape(2) ! shape of the wtcft data
    integer :: fertcft_shape(2) ! shape of the fertcft data
    character(len=*), parameter :: subname = 'dyncrop_init'
    !-----------------------------------------------------------------------
    SHR_ASSERT(bounds%level == BOUNDS_LEVEL_PROC, subname // ': argument must be PROC-level bounds')

    if (masterproc) then
       write(iulog,*) 'Attempting to read crop dynamic landuse data .....'
    end if

    ! Get the year from the START of the timestep; this way, we'll update crop areas
    ! starting after the year boundary. This is consistent with the timing of glacier
    ! updates, and will likely be consistent with the timing of crop updates determined
    ! prognostically, if crop areas are ever determined prognostically rather than
    ! prescribed ahead of time.
    dyncrop_file = dyn_file_type(dyncrop_filename, YEAR_POSITION_START_OF_TIMESTEP)
    call check_dim(dyncrop_file, 'cft', cft_size)
    ! read data PCT_CROP and PCT_CFT corresponding to correct year
    !
    ! Note: if you want to change transient crops so that they are interpolated, rather
    ! than jumping to each year's value on Jan 1 of that year, simply change wtcrop and
    ! wtcft to be of type dyn_var_time_interp_type (rather than
    ! dyn_var_time_uninterp_type), and change the following constructors to construct
    ! variables of dyn_var_time_interp_type. That's all you need to do.
    num_points = (bounds%endg - bounds%begg + 1)
    wtcrop = dyn_var_time_uninterp_type( &
         dyn_file = dyncrop_file, varname=crop_varname, &
         dim1name=grlnd, conversion_factor=100._r8, &
         do_check_sums_equal_1=.false., data_shape=[num_points])
    wtcft_shape = [num_points, cft_size]
    wtcft = dyn_var_time_uninterp_type( &
         dyn_file = dyncrop_file, varname=cft_varname, &
         dim1name=grlnd, conversion_factor=100._r8, &
         do_check_sums_equal_1=.true., data_shape=wtcft_shape)
    fertcft_shape = [num_points, cft_size]
    fertcft = dyn_var_time_uninterp_type( &
         dyn_file = dyncrop_file, varname=fert_varname, &
         dim1name=grlnd, conversion_factor=1._r8, &
         do_check_sums_equal_1=.false., data_shape=fertcft_shape, &
         allow_nodata=.true.)

  end subroutine dyncrop_init

  !-----------------------------------------------------------------------
  subroutine dyncrop_interp(bounds,crop_inst)
    !
    ! !DESCRIPTION:
    ! Get crop cover for model time, when needed.
    !
    ! Sets col_pp%wtlunit and lun%wtgcell for crop landunits.
    !
    ! Note that crop cover currently jumps to its new value at the start of the year.
    ! However, as mentioned above, this behavior can be changed to time interpolation
    ! simply by making wtcrop and wtcft dyn_var_time_interp_type variables rather than
    ! dyn_var_time_uninterp_type.
    !
    ! !USES:
    use CropType          , only : crop_type
    use landunit_varcon   , only : istcrop
    use elm_varpar        , only : cft_lb, cft_ub
    use surfrdUtilsMod    , only : collapse_crop_types
    use subgridWeightsMod , only : set_landunit_weight
    use subgridWeightsMod , only : get_landunit_weight
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds  ! proc-level bounds
    type(crop_type), intent(in) :: crop_inst  ! crop instance for updating annual fertilizer
    !
    ! !LOCAL VARIABLES:
    integer               :: m,p,c,l,g      ! indices
    real(r8), allocatable :: wtcrop_cur(:)  ! current weight of the crop landunit
    real(r8), allocatable :: wtcft_cur(:,:) ! current cft weights
    real(r8), allocatable :: fertcft_cur(:,:) ! current cft fertilizer
    logical , allocatable :: col_set(:)     ! whether we have set the weight for each column
    character(len=*), parameter :: subname = 'dyncrop_interp'
    !-----------------------------------------------------------------------

    SHR_ASSERT(bounds%level == BOUNDS_LEVEL_PROC, subname // ': argument must be PROC-level bounds')

    call dyncrop_file%time_info%set_current_year()

    ! Set new landunit area
    allocate(wtcrop_cur(bounds%begg:bounds%endg))
    call wtcrop%get_current_data(wtcrop_cur)
    do g = bounds%begg, bounds%endg
       call set_landunit_weight(g, istcrop, wtcrop_cur(g))
    end do
    deallocate(wtcrop_cur)

    ! Set new CFT weights
    !
    ! Assumes that memory has been allocated for all CFTs on the crop landunit, and that
    ! each crop is on its own column.
    allocate(wtcft_cur(bounds%begg:bounds%endg, cft_lb:cft_ub))
    call wtcft%get_current_data(wtcft_cur)

    allocate(fertcft_cur(bounds%begg:bounds%endg, cft_lb:cft_ub))
    call fertcft%get_current_data(fertcft_cur)

    call collapse_crop_types(wtcft_cur, fertcft_cur, bounds%begg, bounds%endg, verbose = .false.)

    allocate(col_set(bounds%begc:bounds%endc))
    col_set(:) = .false.

    do p = bounds%begp, bounds%endp
       g = veg_pp%gridcell(p)
       l = veg_pp%landunit(p)
       c = veg_pp%column(p)

       if (lun_pp%itype(l) == istcrop) then
          m = veg_pp%itype(p)

          ! The following assumes there is a single CFT on each crop column. The
          ! error-check with col_set helps ensure this is the case.
          if (col_set(c)) then
             write(iulog,*) subname//' ERROR: attempt to set a column that has already been set.'
             write(iulog,*) 'This may happen if there are multiple crops on a single column.'
             call endrun(decomp_index=c, elmlevel=namec, msg=errMsg(sourcefile, __LINE__))
          end if
          col_pp%wtlunit(c) = wtcft_cur(g,m)
          if (use_crop) then
            crop_inst%fertnitro_patch(p) = fertcft_cur(g,m)
          end if
          col_set(c) = .true.
       end if
    end do

    deallocate(wtcft_cur)
    deallocate(fertcft_cur)
    deallocate(col_set)

  end subroutine dyncrop_interp

end module dyncropFileMod
