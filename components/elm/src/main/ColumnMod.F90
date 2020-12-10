module ColumnMod

  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use elm_varpar     , only : nlevsno, nlevgrnd, nlevlak, nlevslp
  use elm_varcon     , only : spval, ispval
  use ColumnType     , only : col_pp
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

    idx = idx + 1;                                         values(idx) = real( col_pp%landunit (c))
    idx = idx + 1;                                         values(idx) =       col_pp%wtlunit  (c)
    idx = idx + 1;                                         values(idx) = real( col_pp%gridcell (c))
    idx = idx + 1;                                         values(idx) =       col_pp%wtgcell  (c)
    idx = idx + 1;                                         values(idx) = real( col_pp%pfti     (c))
    idx = idx + 1;                                         values(idx) = real( col_pp%pftf     (c))
    idx = idx + 1;                                         values(idx) = real( col_pp%npfts    (c))
    idx = idx + 1;                                         values(idx) = real( col_pp%itype    (c))

    idx = idx + 1; if (col_pp%active(c))                     values(idx) = 1._r8

    idx = idx + 1; if (.not. isnan(col_pp%glc_topo    (c)))  values(idx) = col_pp%glc_topo    (c)
    idx = idx + 1; if (.not. isnan(col_pp%micro_sigma (c)))  values(idx) = col_pp%micro_sigma (c)
    idx = idx + 1; if (.not. isnan(col_pp%n_melt      (c)))  values(idx) = col_pp%n_melt      (c)
    idx = idx + 1; if (.not. isnan(col_pp%topo_slope  (c)))  values(idx) = col_pp%topo_slope  (c)
    idx = idx + 1; if (.not. isnan(col_pp%topo_std    (c)))  values(idx) = col_pp%topo_std    (c)

    idx = idx + 1;                                         values(idx) = real(col_pp%snl   (c))
    idx = idx + 1; if (.not. isnan(col_pp%lakedepth   (c)))  values(idx) =      col_pp%lakedepth(c)

    do j = -nlevsno+1,nlevgrnd
       idx = idx + 1; if (.not. isnan(col_pp%dz(c,j))) values(idx) = col_pp%dz(c,j)
    enddo

    do j = -nlevsno+1,nlevgrnd
       idx = idx + 1; if (.not. isnan(col_pp%z(c,j))) values(idx) = col_pp%z(c,j)
    enddo

    do j = -nlevsno,nlevgrnd
       idx = idx + 1; if (.not. isnan(col_pp%zi(c,j))) values(idx) = col_pp%zi(c,j)
    enddo

    do j = 1,nlevlak
       idx = idx + 1; if (.not. isnan(col_pp%dz_lake(c,j))) values(idx) = col_pp%dz_lake(c,j)
    enddo

    do j = 1,nlevlak
       idx = idx + 1; if (.not. isnan(col_pp%z_lake(c,j))) values(idx) = col_pp%z_lake(c,j)
    enddo

    do j = 1,nlevslp
       idx = idx + 1; if (.not. isnan(col_pp%hslp_p10(c,j))) values(idx) = col_pp%hslp_p10(c,j) 
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
    call NumValuesPerColumn(nvalues)
    if (size(values) /= nvalues) then
       call endrun(msg="ERROR Size of data is incorrect "//errmsg(__FILE__, __LINE__))
    endif

    idx = 0

    idx = idx + 1;                           col_pp%landunit    (c) = int(values(idx))
    idx = idx + 1;                           col_pp%wtlunit     (c) =     values(idx)
    idx = idx + 1;                           col_pp%gridcell    (c) = int(values(idx))
    idx = idx + 1;                           col_pp%wtgcell     (c) = values(idx)
    idx = idx + 1;                           col_pp%pfti        (c) = int(values(idx))
    idx = idx + 1;                           col_pp%pftf        (c) = int(values(idx))
    idx = idx + 1;                           col_pp%npfts       (c) = int(values(idx))
    idx = idx + 1;                           col_pp%itype       (c) = int(values(idx))

    idx = idx + 1; if (values(idx) == 1._r8) col_pp%active      (c) = .true.

    idx = idx + 1;                           col_pp%glc_topo    (c) = values(idx)
    idx = idx + 1;                           col_pp%micro_sigma (c) = values(idx)
    idx = idx + 1;                           col_pp%n_melt      (c) = values(idx)
    idx = idx + 1;                           col_pp%topo_slope  (c) = values(idx)
    idx = idx + 1;                           col_pp%topo_std    (c) = values(idx)

    idx = idx + 1;                           col_pp%snl         (c) = int(values(idx))
    idx = idx + 1;                           col_pp%lakedepth   (c) = values(idx)

    do j = -nlevsno+1,nlevgrnd
       idx = idx + 1; col_pp%dz(c,j) = values(idx)
    enddo

    do j = -nlevsno+1,nlevgrnd
       idx = idx + 1; col_pp%z(c,j) = values(idx)
    enddo

    do j = -nlevsno,nlevgrnd
       idx = idx + 1; col_pp%zi(c,j) = values(idx)
    enddo

    do j = 1,nlevlak
       idx = idx + 1; col_pp%dz_lake(c,j) = values(idx)
    enddo

    do j = 1,nlevlak
       idx = idx + 1; col_pp%z_lake(c,j) = values(idx)
    enddo

    do j = 1,nlevslp
       idx = idx + 1; col_pp%hslp_p10(c,j) = values(idx)
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
