module shr_vcoarsen_mod
  !-------------------------------------------------------------------------------------------
  !
  ! Shared vertical coarsening utilities for E3SM components.
  !
  ! Provides coordinate-agnostic routines for:
  !   1. Overlap-weighted averaging onto coarser vertical layers
  !   2. Level selection by index or nearest coordinate value
  !   3. Category aggregation (sum and weighted average)
  !
  ! All routines are pure computation with no I/O, MPI, or component-specific
  ! dependencies. Each component provides a thin wrapper that builds the
  ! coordinate arrays from its own data structures and calls these routines.
  !
  ! The overlap-weighted averaging algorithm requires monotonically INCREASING
  ! coordinates: coord_iface(k) < coord_iface(k+1) and bounds(j) < bounds(j+1).
  ! This covers pressure (Pa, increasing downward from TOA) and depth (m,
  ! increasing downward from surface).  The overlap formula
  ! max(0, min(c_hi,b_hi) - max(c_lo,b_lo)) assumes c_lo < c_hi and b_lo < b_hi.
  !
  !-------------------------------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none
  private
  save

  ! Overlap-weighted vertical averaging
  public :: shr_vcoarsen_avg        ! single column
  public :: shr_vcoarsen_avg_cols   ! multi-column batch

  ! Level selection
  public :: shr_vcoarsen_select_index    ! extract at a specific level index
  public :: shr_vcoarsen_select_nearest  ! extract at nearest coordinate value

  ! Category aggregation
  public :: shr_vcoarsen_cat_sum     ! sum over categories
  public :: shr_vcoarsen_cat_wgtavg  ! weighted average over categories

contains

  !============================================================================
  subroutine shr_vcoarsen_avg(field, coord_iface, nlev, bounds, n_out, &
                              fillval, field_out)
    !--------------------------------------------------------------------------
    ! Overlap-weighted averaging for a single column.
    !
    ! Computes the weighted average of field(1:nlev) onto n_out target layers
    ! defined by bounds(1:n_out+1). The weight for each input level is the
    ! overlap between the input level extent [coord_iface(k), coord_iface(k+1)]
    ! and the target layer extent [bounds(k_out), bounds(k_out+1)].
    !
    ! Works for any monotonic coordinate. Caller provides:
    !   - coord_iface: interface values of the input levels (size nlev+1)
    !   - bounds: boundaries of the target coarsened layers (size n_out+1)
    !   - fillval: value to use when no input levels overlap a target layer
    !
    ! For atmosphere: coord_iface = pressure interfaces (Pa), bounds = pressure bounds
    ! For ocean: coord_iface = depth interfaces (m), bounds = depth bounds
    ! For land: coord_iface = soil depth interfaces (m), bounds = depth bounds
    !--------------------------------------------------------------------------
    integer,  intent(in)  :: nlev              ! number of input vertical levels
    integer,  intent(in)  :: n_out             ! number of output layers
    real(r8), intent(in)  :: field(nlev)       ! input field values
    real(r8), intent(in)  :: coord_iface(nlev+1) ! interface coordinates
    real(r8), intent(in)  :: bounds(n_out+1)   ! target layer boundaries
    real(r8), intent(in)  :: fillval           ! fill value for empty layers
    real(r8), intent(out) :: field_out(n_out)  ! output coarsened values

    real(r8) :: b_lo, b_hi, c_lo, c_hi, overlap
    real(r8) :: numerator, denominator
    integer  :: k, kout

    do kout = 1, n_out
      b_lo = bounds(kout)
      b_hi = bounds(kout + 1)

      numerator   = 0.0_r8
      denominator = 0.0_r8

      do k = 1, nlev
        c_lo = coord_iface(k)
        c_hi = coord_iface(k + 1)

        ! Overlap between [c_lo,c_hi] and [b_lo,b_hi]
        overlap = max(0.0_r8, min(c_hi, b_hi) - max(c_lo, b_lo))

        if (overlap > 0.0_r8) then
          ! Skip NaN values (NaN /= NaN) and fill values
          if (field(k) /= field(k) .or. field(k) == fillval) cycle
          numerator   = numerator   + field(k) * overlap
          denominator = denominator + overlap
        end if
      end do

      if (denominator > 0.0_r8) then
        field_out(kout) = numerator / denominator
      else
        field_out(kout) = fillval
      end if
    end do

  end subroutine shr_vcoarsen_avg

  !============================================================================
  subroutine shr_vcoarsen_avg_cols(field, coord_iface, ncol, nlev, bounds, &
                                    n_out, fillval, field_out)
    !--------------------------------------------------------------------------
    ! Overlap-weighted averaging for multiple columns (batch version).
    !
    ! Same algorithm as shr_vcoarsen_avg applied to each of ncol columns.
    ! Suitable for components with uniform level counts across columns
    ! (e.g., EAM atmosphere, ELM soil).
    !
    ! coord_iface can vary per column (e.g., pressure interfaces in EAM)
    ! or be replicated (e.g., soil interfaces in ELM are the same for all columns).
    !--------------------------------------------------------------------------
    integer,  intent(in)  :: ncol                      ! number of columns
    integer,  intent(in)  :: nlev                      ! number of input levels
    integer,  intent(in)  :: n_out                     ! number of output layers
    real(r8), intent(in)  :: field(ncol, nlev)         ! input field values
    real(r8), intent(in)  :: coord_iface(ncol, nlev+1) ! interface coordinates
    real(r8), intent(in)  :: bounds(n_out+1)           ! target layer boundaries
    real(r8), intent(in)  :: fillval                   ! fill value for empty layers
    real(r8), intent(out) :: field_out(ncol, n_out)    ! output coarsened values

    real(r8) :: b_lo, b_hi, overlap
    real(r8) :: numerator(ncol), denominator(ncol)
    integer  :: i, k, kout

    do kout = 1, n_out
      b_lo = bounds(kout)
      b_hi = bounds(kout + 1)

      numerator(:)   = 0.0_r8
      denominator(:) = 0.0_r8

      ! k outer, i inner: field(i,k) with varying i is contiguous in memory
      do k = 1, nlev
        do i = 1, ncol
          overlap = max(0.0_r8, min(coord_iface(i, k+1), b_hi) - &
                                max(coord_iface(i, k),   b_lo))

          if (overlap > 0.0_r8) then
            ! Skip NaN values (NaN /= NaN) and fill values
            if (field(i, k) /= field(i, k) .or. field(i, k) == fillval) cycle
            numerator(i)   = numerator(i)   + field(i, k) * overlap
            denominator(i) = denominator(i) + overlap
          end if
        end do
      end do

      do i = 1, ncol
        if (denominator(i) > 0.0_r8) then
          field_out(i, kout) = numerator(i) / denominator(i)
        else
          field_out(i, kout) = fillval
        end if
      end do
    end do

  end subroutine shr_vcoarsen_avg_cols

  !============================================================================
  subroutine shr_vcoarsen_select_index(field, ncol, nlev, k_sel, nlev_max, &
                                        fillval, field_out)
    !--------------------------------------------------------------------------
    ! Extract the field value at a specific vertical level index.
    !
    ! For each column i, returns field(i, k_sel) if k_sel <= nlev_max(i),
    ! otherwise returns fillval.
    !
    ! Output naming convention: base_at_L{k_sel} (e.g., "U_at_L5")
    !--------------------------------------------------------------------------
    integer,  intent(in)  :: ncol              ! number of columns
    integer,  intent(in)  :: nlev              ! number of input levels
    real(r8), intent(in)  :: field(ncol, nlev) ! input field values
    integer,  intent(in)  :: k_sel             ! level index to select
    integer,  intent(in)  :: nlev_max(ncol)    ! max valid level per column
    real(r8), intent(in)  :: fillval           ! fill value for invalid columns
    real(r8), intent(out) :: field_out(ncol)   ! output selected values

    integer :: i

    if (k_sel < 1 .or. k_sel > nlev) then
      field_out(1:ncol) = fillval
      return
    end if

    do i = 1, ncol
      if (k_sel <= nlev_max(i)) then
        field_out(i) = field(i, k_sel)
      else
        field_out(i) = fillval
      end if
    end do

  end subroutine shr_vcoarsen_select_index

  !============================================================================
  subroutine shr_vcoarsen_select_nearest(field, coord_mid, ncol, nlev, &
                                          target_val, nlev_max, fillval, &
                                          field_out)
    !--------------------------------------------------------------------------
    ! Linearly interpolate the field to a target coordinate value.
    !
    ! For each column, finds the two adjacent levels that bracket target_val
    ! and linearly interpolates between them. If the target is outside the
    ! coordinate range, the nearest boundary level value is returned (no
    ! extrapolation). Works for any monotonic coordinate direction.
    !
    ! If one bracketing level has a fill/NaN value, the other valid level
    ! is used. If both are invalid, fillval is returned.
    !
    ! Output naming convention: base_at_P{target} (e.g., "U_at_P850")
    ! The P prefix is component-specific (P for pressure, D for depth, etc.)
    !--------------------------------------------------------------------------
    integer,  intent(in)  :: ncol                  ! number of columns
    integer,  intent(in)  :: nlev                  ! number of input levels
    real(r8), intent(in)  :: field(ncol, nlev)     ! input field values
    real(r8), intent(in)  :: coord_mid(ncol, nlev) ! level midpoint coordinates
    real(r8), intent(in)  :: target_val            ! target coordinate value
    integer,  intent(in)  :: nlev_max(ncol)        ! max valid level per column
    real(r8), intent(in)  :: fillval               ! fill value for invalid columns
    real(r8), intent(out) :: field_out(ncol)       ! output selected values

    integer  :: i, k, k_lo
    real(r8) :: d_total, w, f_lo, f_hi
    logical  :: lo_valid, hi_valid

    do i = 1, ncol
      if (nlev_max(i) < 1) then
        field_out(i) = fillval
        cycle
      end if

      ! Single valid level: no interpolation possible
      if (nlev_max(i) == 1) then
        f_lo = field(i, 1)
        if (f_lo /= f_lo .or. f_lo == fillval) then
          field_out(i) = fillval
        else
          field_out(i) = f_lo
        end if
        cycle
      end if

      ! Find bracket: k_lo such that target is between coord(k_lo) and coord(k_lo+1).
      ! The sign-product test works for both increasing and decreasing coordinates.
      k_lo = 0
      do k = 1, nlev_max(i) - 1
        if ((coord_mid(i, k) - target_val) * &
            (coord_mid(i, k+1) - target_val) <= 0.0_r8) then
          k_lo = k
          exit
        end if
      end do

      if (k_lo < 1) then
        ! Target outside coordinate range — use nearest boundary (no extrapolation)
        if (abs(coord_mid(i, 1) - target_val) <= &
            abs(coord_mid(i, nlev_max(i)) - target_val)) then
          f_lo = field(i, 1)
        else
          f_lo = field(i, nlev_max(i))
        end if
        if (f_lo /= f_lo .or. f_lo == fillval) then
          field_out(i) = fillval
        else
          field_out(i) = f_lo
        end if
      else
        ! Linear interpolation between k_lo and k_lo+1
        f_lo = field(i, k_lo)
        f_hi = field(i, k_lo + 1)
        lo_valid = (f_lo == f_lo .and. f_lo /= fillval)
        hi_valid = (f_hi == f_hi .and. f_hi /= fillval)

        if (lo_valid .and. hi_valid) then
          d_total = coord_mid(i, k_lo+1) - coord_mid(i, k_lo)
          if (abs(d_total) > 0.0_r8) then
            w = (target_val - coord_mid(i, k_lo)) / d_total
            field_out(i) = f_lo * (1.0_r8 - w) + f_hi * w
          else
            field_out(i) = f_lo
          end if
        else if (lo_valid) then
          field_out(i) = f_lo
        else if (hi_valid) then
          field_out(i) = f_hi
        else
          field_out(i) = fillval
        end if
      end if
    end do

  end subroutine shr_vcoarsen_select_nearest

  !============================================================================
  subroutine shr_vcoarsen_cat_sum(field, ncol, ncat, field_out)
    !--------------------------------------------------------------------------
    ! Sum field values over categories.
    !
    ! For each column: field_out(i) = sum(field(i, 1:ncat))
    !
    ! Applicable to sea ice category aggregation (e.g., total ice area,
    ! total ice volume, total snow volume).
    !--------------------------------------------------------------------------
    integer,  intent(in)  :: ncol              ! number of columns/cells
    integer,  intent(in)  :: ncat              ! number of categories
    real(r8), intent(in)  :: field(ncol, ncat) ! input field values per category
    real(r8), intent(out) :: field_out(ncol)   ! output aggregated values

    integer :: i, k

    do i = 1, ncol
      field_out(i) = 0.0_r8
      do k = 1, ncat
        field_out(i) = field_out(i) + field(i, k)
      end do
    end do

  end subroutine shr_vcoarsen_cat_sum

  !============================================================================
  subroutine shr_vcoarsen_cat_wgtavg(field, weights, ncol, ncat, fillval, &
                                      field_out)
    !--------------------------------------------------------------------------
    ! Weighted average of field values over categories.
    !
    ! For each column:
    !   field_out(i) = sum(field(i,k) * weights(i,k)) / sum(weights(i,k))
    !
    ! If total weight is zero, field_out(i) = fillval.
    !
    ! Applicable to sea ice area-weighted averages (e.g., mean surface
    ! temperature weighted by ice area fraction per category).
    !--------------------------------------------------------------------------
    integer,  intent(in)  :: ncol                  ! number of columns/cells
    integer,  intent(in)  :: ncat                  ! number of categories
    real(r8), intent(in)  :: field(ncol, ncat)     ! input field values
    real(r8), intent(in)  :: weights(ncol, ncat)   ! weights per category
    real(r8), intent(in)  :: fillval               ! fill value when total weight is zero
    real(r8), intent(out) :: field_out(ncol)       ! output weighted averages

    real(r8) :: numerator, denominator
    integer  :: i, k

    do i = 1, ncol
      numerator   = 0.0_r8
      denominator = 0.0_r8

      do k = 1, ncat
        if (weights(i, k) > 0.0_r8) then
          numerator   = numerator   + field(i, k) * weights(i, k)
          denominator = denominator + weights(i, k)
        end if
      end do

      if (denominator > 0.0_r8) then
        field_out(i) = numerator / denominator
      else
        field_out(i) = fillval
      end if
    end do

  end subroutine shr_vcoarsen_cat_wgtavg

end module shr_vcoarsen_mod
