module ColumnMod

  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use clm_varpar     , only : nlevsno, nlevgrnd, nlevlak
  use clm_varcon     , only : spval, ispval
  use ColumnType     , only : column_type, col
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  public :: GetValuesForColumn
  public :: SetValuesForColumn
  public :: NumValuesPerColumn

  !-----------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine GetValuesForColumn(c, values)
    !
    ! !DESCRIPTION:
    ! Returns all values for the c-th column
    !
    ! NOTE: The order of the values retrieved/stored in
    !       GetValuesForColumn/SetValuesForColumn should be same.
    !
    ! !USES
    use abortutils       , only : endrun
    use shr_infnan_mod   , only : isnan => shr_infnan_isnan
    use shr_log_mod      , only : errMsg => shr_log_errMsg
    !
    ! !ARGUMENTS:
    integer, intent(in)            :: c
    real(r8), pointer, intent(out) :: values(:)
    !
    ! !LOCAL VARIABLES:
    integer :: idx
    integer :: j
    integer :: nvalues

    ! Check the size of data
    call NumValuesPerColumn(nvalues)
    if (size(values) /= nvalues) then
       call endrun(msg="ERROR Size of data is incorrect "//errmsg(__FILE__, __LINE__))
    endif

    idx = 0

    idx = idx + 1;                                         values(idx) = real( col%landunit (c))
    idx = idx + 1;                                         values(idx) =       col%wtlunit  (c)
    idx = idx + 1;                                         values(idx) = real( col%gridcell (c))
    idx = idx + 1;                                         values(idx) =       col%wtgcell  (c)
    idx = idx + 1;                                         values(idx) = real( col%pfti     (c))
    idx = idx + 1;                                         values(idx) = real( col%pftf     (c))
    idx = idx + 1;                                         values(idx) = real( col%npfts    (c))
    idx = idx + 1;                                         values(idx) = real( col%itype    (c))

    idx = idx + 1; if (col%active(c))                     values(idx) = 1._r8

    idx = idx + 1; if (.not. isnan(col%glc_topo    (c)))  values(idx) = col%glc_topo    (c)
    idx = idx + 1; if (.not. isnan(col%micro_sigma (c)))  values(idx) = col%micro_sigma (c)
    idx = idx + 1; if (.not. isnan(col%n_melt      (c)))  values(idx) = col%n_melt      (c)
    idx = idx + 1; if (.not. isnan(col%topo_slope  (c)))  values(idx) = col%topo_slope  (c)
    idx = idx + 1; if (.not. isnan(col%topo_std    (c)))  values(idx) = col%topo_std    (c)

    idx = idx + 1;                                         values(idx) = real(col%snl   (c))
    idx = idx + 1; if (.not. isnan(col%lakedepth   (c)))  values(idx) =      col%lakedepth(c)

    do j = -nlevsno+1,nlevgrnd
       idx = idx + 1; if (.not. isnan(col%dz(c,j))) values(idx) = col%dz(c,j)
    enddo

    do j = -nlevsno+1,nlevgrnd
       idx = idx + 1; if (.not. isnan(col%z(c,j))) values(idx) = col%z(c,j)
    enddo

    do j = -nlevsno,nlevgrnd
       idx = idx + 1; if (.not. isnan(col%zi(c,j))) values(idx) = col%zi(c,j)
    enddo

    do j = 1,nlevlak
       idx = idx + 1; if (.not. isnan(col%dz_lake(c,j))) values(idx) = col%dz_lake(c,j)
    enddo

    do j = 1,nlevlak
       idx = idx + 1; if (.not. isnan(col%z_lake(c,j))) values(idx) = col%z_lake(c,j)
    enddo

  end subroutine GetValuesForColumn

  !------------------------------------------------------------------------
  subroutine SetValuesForColumn(c, values)
    !
    ! !DESCRIPTION:
    ! Sets all values for the c-th column.
    !
    ! NOTE: The order of the values retrieved/stored in
    !       GetValuesForColumn/SetValuesForColumn should be same.
    !
    ! !USES
    use abortutils       , only : endrun
    use shr_infnan_mod   , only : isnan => shr_infnan_isnan
    use shr_log_mod      , only : errMsg => shr_log_errMsg
    !
    ! !ARGUMENTS:
    integer, intent(in)            :: c
    real(r8), pointer, intent(in)  :: values(:)
    !
    ! !LOCAL VARIABLES:
    integer :: idx
    integer :: j
    integer :: nvalues

    ! Check the size of data
    if (size(values) /= nvalues) then
       call endrun(msg="ERROR Size of data is incorrect "//errmsg(__FILE__, __LINE__))
    endif

    idx = 0

    idx = idx + 1;                           col%landunit    (c) = int(values(idx))
    idx = idx + 1;                           col%wtlunit     (c) =     values(idx)
    idx = idx + 1;                           col%gridcell    (c) = int(values(idx))
    idx = idx + 1;                           col%wtgcell     (c) = values(idx)
    idx = idx + 1;                           col%pfti        (c) = int(values(idx))
    idx = idx + 1;                           col%pftf        (c) = int(values(idx))
    idx = idx + 1;                           col%npfts       (c) = int(values(idx))
    idx = idx + 1;                           col%itype       (c) = int(values(idx))

    idx = idx + 1; if (values(idx) == 1._r8) col%active      (c) = .true.

    idx = idx + 1;                           col%glc_topo    (c) = values(idx)
    idx = idx + 1;                           col%micro_sigma (c) = values(idx)
    idx = idx + 1;                           col%n_melt      (c) = values(idx)
    idx = idx + 1;                           col%topo_slope  (c) = values(idx)
    idx = idx + 1;                           col%topo_std    (c) = values(idx)

    idx = idx + 1;                           col%snl         (c) = int(values(idx))
    idx = idx + 1;                           col%lakedepth   (c) = values(idx)

    do j = -nlevsno+1,nlevgrnd
       idx = idx + 1; col%dz(c,j) = values(idx)
    enddo

    do j = -nlevsno+1,nlevgrnd
       idx = idx + 1; col%z(c,j) = values(idx)
    enddo

    do j = -nlevsno,nlevgrnd
       idx = idx + 1; col%zi(c,j) = values(idx)
    enddo

    do j = 1,nlevlak
       idx = idx + 1; col%dz_lake(c,j) = values(idx)
    enddo

    do j = 1,nlevlak
       idx = idx + 1; col%z_lake(c,j) = values(idx)
    enddo

  end subroutine SetValuesForColumn

  !------------------------------------------------------------------------
  subroutine NumValuesPerColumn(nvalues)
    !
    ! !DESCRIPTION:
    ! Returns the number of values for each column within the column_type
    !
    ! !ARGUMENTS:
    integer, intent(out)     :: nvalues

    nvalues = 17 + 2*(nlevgrnd+nlevsno+1) + (nlevgrnd+nlevsno) + 2*nlevlak

  end subroutine NumValuesPerColumn


end module ColumnMod
