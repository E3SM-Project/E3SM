module histFieldsMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing subroutines to initialize history fields for each
  ! data type.
  !
  ! !USES:
  use clm_varcon     , only : spval, ispval
  !use histFileMod    , only : hist_addfld1d, hist_addfld2d, no_snow_normal
  use ColumnType     , only : col_es
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: col_es_init_hist        ! Initialize history fields for the column energy state data type
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine col_es_init_hist(begc, endc)
    !
    ! !DESCRIPTION:
    ! Initialize history fields for the column energy state data type
    !
    ! !ARGUMENTS:
    integer, intent(in)                   :: begc, endc
    !-----------------------------------------------------------------------

    !col_es%t_h2osfc(begc:endc) = spval
    !call hist_addfld1d (fname='TH2OSFC',  units='K',  &
    !     avgflag='A', long_name='surface water temperature', &
    !     ptr_col=col_es%t_h2osfc)
         
  end subroutine col_es_init_hist

end module histFieldsMod