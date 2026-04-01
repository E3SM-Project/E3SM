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
  ! The overlap-weighted averaging algorithm works for any monotonic vertical
  ! coordinate (pressure, depth, etc.) as long as coord_iface and bounds use
  ! the same convention (both increasing or both decreasing).
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
    real(r8), intent(in)  :: field(nlev)       ! input field values
    real(r8), intent(in)  :: coord_iface(nlev+1) ! interface coordinates
    real(r8), intent(in)  :: bounds(n_out+1)   ! target layer boundaries
    integer,  intent(in)  :: n_out             ! number of output layers
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
    real(r8), intent(in)  :: field(ncol, nlev)         ! input field values
    real(r8), intent(in)  :: coord_iface(ncol, nlev+1) ! interface coordinates
    real(r8), intent(in)  :: bounds(n_out+1)           ! target layer boundaries
    integer,  intent(in)  :: n_out                     ! number of output layers
    real(r8), intent(in)  :: fillval                   ! fill value for empty layers
    real(r8), intent(out) :: field_out(ncol, n_out)    ! output coarsened values

    real(r8) :: b_lo, b_hi, c_lo, c_hi, overlap
    real(r8) :: numerator, denominator
    integer  :: i, k, kout

    do kout = 1, n_out
      b_lo = bounds(kout)
      b_hi = bounds(kout + 1)

      do i = 1, ncol
        numerator   = 0.0_r8
        denominator = 0.0_r8

        do k = 1, nlev
          c_lo = coord_iface(i, k)
          c_hi = coord_iface(i, k + 1)

          overlap = max(0.0_r8, min(c_hi, b_hi) - max(c_lo, b_lo))

          if (overlap > 0.0_r8) then
            numerator   = numerator   + field(i, k) * overlap
            denominator = denominator + overlap
          end if
        end do

        if (denominator > 0.0_r8) then
          field_out(i, kout) = numerator / denominator
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
    ! Extract the field value at the level nearest to a target coordinate value.
    !
    ! For each column, finds the level k (within 1..nlev_max(i)) where
    ! |coord_mid(i,k) - target_val| is minimized, and returns field(i,k).
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

    integer  :: i, k, k_best
    real(r8) :: dist, dist_best

    do i = 1, ncol
      if (nlev_max(i) < 1) then
        field_out(i) = fillval
        cycle
      end if

      k_best = 1
      dist_best = abs(coord_mid(i, 1) - target_val)

      do k = 2, nlev_max(i)
        dist = abs(coord_mid(i, k) - target_val)
        if (dist < dist_best) then
          dist_best = dist
          k_best = k
        end if
      end do

      field_out(i) = field(i, k_best)
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
