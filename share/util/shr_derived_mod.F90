module shr_derived_mod
  !-------------------------------------------------------------------------------------------
  !
  ! Shared derived field utilities for E3SM components.
  !
  ! Provides component-independent routines for:
  !   1. Parsing derived field definitions (e.g., "TOTAL_WATER=Q+CLDICE+CLDLIQ+QRAIN")
  !   2. Evaluating parsed expressions given field data arrays
  !   3. Computing tendencies from cached previous-step values
  !
  ! Supports:
  !   - Binary operators: +, -, *, / between fields and/or numeric constants
  !   - Strict left-to-right evaluation (no operator precedence)
  !   - Auto-detection of numeric constants vs field names
  !
  ! All routines are pure computation with no I/O, MPI, or component-specific
  ! dependencies. Each component provides a thin wrapper that handles its own
  ! field lookup, namelist I/O, and history system integration.
  !
  ! Expression syntax:
  !   OUTPUT_NAME=OPERAND1 op OPERAND2 op OPERAND3 ...
  !   where op is +, -, *, or /
  !   and each operand is either a field name or a numeric constant.
  !   Numeric constants are auto-detected (e.g., 9.8, 1000.0, 8.64e7).
  !
  ! Examples:
  !   TOTAL_WATER=Q+CLDICE+CLDLIQ+QRAIN    (sum of 4 fields)
  !   HGTsfc=PHIS/9.80616                    (field / constant)
  !   PRECT_mm=PRECT*8.64e7                  (field * constant)
  !   TOTAL_WATER_g=TOTAL_WATER*1000.0       (chained derived * constant)
  !
  !-------------------------------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none
  private
  save

  ! Public parameters
  integer, parameter, public :: shr_derived_max_operands = 20
  integer, parameter, public :: shr_derived_max_namelen  = 34
  integer, parameter, public :: shr_derived_max_deflen   = 256

  ! Public types
  public :: shr_derived_operand_t
  public :: shr_derived_expr_t

  ! Public routines
  public :: shr_derived_parse
  public :: shr_derived_eval
  public :: shr_derived_tend
  public :: shr_derived_is_number

  type :: shr_derived_operand_t
    logical                                    :: is_constant = .false.
    real(r8)                                   :: const_val   = 0.0_r8
    character(len=shr_derived_max_namelen)      :: field_name  = ''
    character(len=1)                            :: op          = '+'
  end type shr_derived_operand_t

  type :: shr_derived_expr_t
    character(len=shr_derived_max_namelen)      :: output_name = ''
    integer                                     :: n_operands  = 0
    type(shr_derived_operand_t)                 :: operands(shr_derived_max_operands)
    character(len=shr_derived_max_deflen)       :: long_name   = ''
  end type shr_derived_expr_t

contains

  !============================================================================
  subroutine shr_derived_parse(defstr, expr, ierr)
    !--------------------------------------------------------------------------
    ! Parse a derived field definition string.
    !
    ! Format: "OUTPUT_NAME=OPERAND1+OPERAND2-OPERAND3*OPERAND4/OPERAND5"
    !
    ! Each operand is either a field name or a numeric constant.
    ! Constants are auto-detected by attempting a Fortran internal read.
    ! The first operand has an implicit '+' operator.
    ! Operators bind the NEXT operand (e.g., A+B means: A with '+', B with '+').
    ! Evaluation is strict left-to-right with no operator precedence.
    !
    ! ierr = 0 on success, nonzero on error.
    !--------------------------------------------------------------------------
    character(len=*),          intent(in)  :: defstr
    type(shr_derived_expr_t), intent(out) :: expr
    integer,                  intent(out) :: ierr

    integer :: eq_pos, start_pos, i, n, rhs_len
    character(len=shr_derived_max_deflen) :: rhs
    character(len=shr_derived_max_namelen) :: token
    character(len=1) :: ch, next_op
    real(r8) :: val

    ierr = 0
    expr%n_operands = 0
    expr%output_name = ''
    expr%long_name = ''

    ! Find '=' separator
    eq_pos = index(defstr, '=')
    if (eq_pos < 2) then
      ierr = 1
      return
    end if

    ! Extract output name
    expr%output_name = adjustl(defstr(1:eq_pos-1))

    ! Extract RHS
    rhs = adjustl(defstr(eq_pos+1:))
    rhs_len = len_trim(rhs)
    if (rhs_len == 0) then
      ierr = 2
      return
    end if

    ! Parse operator-separated operands
    n = 0
    start_pos = 1
    next_op = '+'  ! first operand has implicit '+'

    do
      ! Find next operator
      i = start_pos
      do while (i <= rhs_len)
        ch = rhs(i:i)
        ! Skip if this could be part of scientific notation (e.g., 1.5e-3)
        if ((ch == '+' .or. ch == '-') .and. i > start_pos) then
          ! Check if preceded by 'e' or 'E' (scientific notation)
          if (rhs(i-1:i-1) == 'e' .or. rhs(i-1:i-1) == 'E') then
            ! Only treat as scientific notation if char before e/E is a digit
            if (i > start_pos + 1) then
              if (rhs(i-2:i-2) >= '0' .and. rhs(i-2:i-2) <= '9') then
                i = i + 1
                cycle
              end if
            end if
          end if
        end if
        if (ch == '+' .or. ch == '-' .or. ch == '*' .or. ch == '/') then
          if (i > start_pos) exit  ! found operator after token
        end if
        i = i + 1
      end do

      ! Extract token
      n = n + 1
      if (n > shr_derived_max_operands) then
        ierr = 3
        return
      end if

      if (i <= rhs_len) then
        ! Operator found at position i
        token = adjustl(rhs(start_pos:i-1))
        expr%operands(n)%op = next_op
        next_op = rhs(i:i)
        start_pos = i + 1
      else
        ! No more operators, rest is last token
        token = adjustl(rhs(start_pos:rhs_len))
        expr%operands(n)%op = next_op
      end if

      ! Determine if token is a constant or field name
      if (shr_derived_is_number(token, val)) then
        expr%operands(n)%is_constant = .true.
        expr%operands(n)%const_val = val
        expr%operands(n)%field_name = ''
      else
        expr%operands(n)%is_constant = .false.
        expr%operands(n)%const_val = 0.0_r8
        expr%operands(n)%field_name = token
      end if

      if (i > rhs_len) exit
    end do

    expr%n_operands = n

    ! Build long_name from the RHS expression
    expr%long_name = trim(rhs)

  end subroutine shr_derived_parse

  !============================================================================
  subroutine shr_derived_eval(expr, fields, ncol, nlev, max_col, result)
    !--------------------------------------------------------------------------
    ! Evaluate a parsed expression given field data arrays.
    !
    ! fields(max_col, nlev, n_operands): caller provides field data for each
    !   non-constant operand. For constant operands, the corresponding slice
    !   is ignored (can be uninitialized).
    !
    ! Evaluation is strict left-to-right:
    !   A+B*C  evaluates as  (A+B)*C
    !   A*B+C  evaluates as  (A*B)+C
    !
    ! Division by zero produces 0.0.
    !--------------------------------------------------------------------------
    type(shr_derived_expr_t), intent(in)  :: expr
    integer,                  intent(in)  :: ncol
    integer,                  intent(in)  :: nlev
    integer,                  intent(in)  :: max_col
    real(r8),                 intent(in)  :: fields(max_col, nlev, expr%n_operands)
    real(r8),                 intent(out) :: result(max_col, nlev)

    real(r8) :: src(max_col, nlev)
    integer  :: n

    result(1:ncol, 1:nlev) = 0.0_r8

    do n = 1, expr%n_operands
      ! Load source: constant or field
      if (expr%operands(n)%is_constant) then
        src(1:ncol, 1:nlev) = expr%operands(n)%const_val
      else
        src(1:ncol, 1:nlev) = fields(1:ncol, 1:nlev, n)
      end if

      ! Apply operator
      select case (expr%operands(n)%op)
      case ('+')
        if (n == 1) then
          result(1:ncol, 1:nlev) = src(1:ncol, 1:nlev)
        else
          result(1:ncol, 1:nlev) = result(1:ncol, 1:nlev) + src(1:ncol, 1:nlev)
        end if
      case ('-')
        if (n == 1) then
          result(1:ncol, 1:nlev) = -src(1:ncol, 1:nlev)
        else
          result(1:ncol, 1:nlev) = result(1:ncol, 1:nlev) - src(1:ncol, 1:nlev)
        end if
      case ('*')
        if (n == 1) then
          result(1:ncol, 1:nlev) = src(1:ncol, 1:nlev)
        else
          result(1:ncol, 1:nlev) = result(1:ncol, 1:nlev) * src(1:ncol, 1:nlev)
        end if
      case ('/')
        if (n == 1) then
          result(1:ncol, 1:nlev) = src(1:ncol, 1:nlev)
        else
          where (src(1:ncol, 1:nlev) /= 0.0_r8)
            result(1:ncol, 1:nlev) = result(1:ncol, 1:nlev) / src(1:ncol, 1:nlev)
          elsewhere
            result(1:ncol, 1:nlev) = 0.0_r8
          end where
        end if
      end select
    end do

  end subroutine shr_derived_eval

  !============================================================================
  subroutine shr_derived_tend(curr, prev, ncol, nlev, max_col, dt, tend)
    !--------------------------------------------------------------------------
    ! Compute tendency as (curr - prev) / dt.
    !
    ! This is the pure-computation part of tendency calculation.
    ! The component wrapper is responsible for caching prev values between
    ! timesteps and managing first-timestep logic.
    !
    ! Works for both 2D (nlev=1) and 3D (nlev=pver) fields.
    !--------------------------------------------------------------------------
    integer,  intent(in)  :: ncol
    integer,  intent(in)  :: nlev
    integer,  intent(in)  :: max_col
    real(r8), intent(in)  :: curr(max_col, nlev)   ! current timestep values
    real(r8), intent(in)  :: prev(max_col, nlev)   ! previous timestep values
    real(r8), intent(in)  :: dt                     ! timestep in seconds
    real(r8), intent(out) :: tend(max_col, nlev)    ! output tendency

    if (dt > 0.0_r8) then
      tend(1:ncol, 1:nlev) = (curr(1:ncol, 1:nlev) - prev(1:ncol, 1:nlev)) / dt
    else
      tend(1:ncol, 1:nlev) = 0.0_r8
    end if

  end subroutine shr_derived_tend

  !============================================================================
  function shr_derived_is_number(token, val) result(is_num)
    !--------------------------------------------------------------------------
    ! Check if a string token represents a numeric constant.
    ! Returns .true. and sets val if the token parses as a real number.
    ! Handles integers, decimals, and scientific notation (e.g., 8.64e7).
    !--------------------------------------------------------------------------
    character(len=*), intent(in)  :: token
    real(r8),         intent(out) :: val
    logical :: is_num

    integer :: ios

    val = 0.0_r8
    is_num = .false.

    if (len_trim(token) == 0) return

    read(token, *, iostat=ios) val
    is_num = (ios == 0)

  end function shr_derived_is_number

end module shr_derived_mod
