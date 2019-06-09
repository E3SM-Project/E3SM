    !-------------------------------------------------------------------
    ! Overloaded definitions for (+):
    !

    ELEMENTAL FUNCTION add_rpe (x) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        TYPE(rpe_var) :: z
        z%sbits = significand_bits(x)
        z = +(x%val)
    END FUNCTION add_rpe

    ELEMENTAL FUNCTION add_rpe_rpe (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        TYPE(rpe_var) :: z
        z%sbits = MAX(significand_bits(x), significand_bits(y))
        z = x%val + y%val
    END FUNCTION add_rpe_rpe

    ELEMENTAL FUNCTION add_rpe_integer (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        INTEGER, INTENT(IN) :: y
        TYPE(rpe_var) :: z
        z%sbits = MAX(significand_bits(x), significand_bits(y))
        z = x%val + y
    END FUNCTION add_rpe_integer

    ELEMENTAL FUNCTION add_rpe_long (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        INTEGER(KIND=8), INTENT(IN) :: y
        TYPE(rpe_var) :: z
        z%sbits = MAX(significand_bits(x), significand_bits(y))
        z = x%val + y
    END FUNCTION add_rpe_long

    ELEMENTAL FUNCTION add_rpe_real (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        REAL(KIND=RPE_REAL_KIND), INTENT(IN) :: y
        TYPE(rpe_var) :: z
        z%sbits = MAX(significand_bits(x), significand_bits(y))
        z = x%val + y
    END FUNCTION add_rpe_real

    ELEMENTAL FUNCTION add_rpe_realalt (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        REAL(KIND=RPE_ALTERNATE_KIND), INTENT(IN) :: y
        TYPE(rpe_var) :: z
        z%sbits = MAX(significand_bits(x), significand_bits(y))
        z = x%val + y
    END FUNCTION add_rpe_realalt

    ELEMENTAL FUNCTION add_integer_rpe (x, y) RESULT (z)
        INTEGER, INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        TYPE(rpe_var) :: z
        z%sbits = MAX(significand_bits(x), significand_bits(y))
        z = x + y%val
    END FUNCTION add_integer_rpe

    ELEMENTAL FUNCTION add_long_rpe (x, y) RESULT (z)
        INTEGER(KIND=8), INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        TYPE(rpe_var) :: z
        z%sbits = MAX(significand_bits(x), significand_bits(y))
        z = x + y%val
    END FUNCTION add_long_rpe

    ELEMENTAL FUNCTION add_real_rpe (x, y) RESULT (z)
        REAL(KIND=RPE_REAL_KIND), INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        TYPE(rpe_var) :: z
        z%sbits = MAX(significand_bits(x), significand_bits(y))
        z = x + y%val
    END FUNCTION add_real_rpe

    ELEMENTAL FUNCTION add_realalt_rpe (x, y) RESULT (z)
        REAL(KIND=RPE_ALTERNATE_KIND), INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        TYPE(rpe_var) :: z
        z%sbits = MAX(significand_bits(x), significand_bits(y))
        z = x + y%val
    END FUNCTION add_realalt_rpe

    !-------------------------------------------------------------------
    ! Overloaded definitions for (-):
    !

    ELEMENTAL FUNCTION sub_rpe (x) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        TYPE(rpe_var) :: z
        z%sbits = significand_bits(x)
        z = -(x%val)
    END FUNCTION sub_rpe

    ELEMENTAL FUNCTION sub_rpe_rpe (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        TYPE(rpe_var) :: z
        z%sbits = MAX(significand_bits(x), significand_bits(y))
        z = x%val - y%val
    END FUNCTION sub_rpe_rpe

    ELEMENTAL FUNCTION sub_rpe_integer (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        INTEGER, INTENT(IN) :: y
        TYPE(rpe_var) :: z
        z%sbits = MAX(significand_bits(x), significand_bits(y))
        z = x%val - y
    END FUNCTION sub_rpe_integer

    ELEMENTAL FUNCTION sub_rpe_long (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        INTEGER(KIND=8), INTENT(IN) :: y
        TYPE(rpe_var) :: z
        z%sbits = MAX(significand_bits(x), significand_bits(y))
        z = x%val - y
    END FUNCTION sub_rpe_long

    ELEMENTAL FUNCTION sub_rpe_real (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        REAL(KIND=RPE_REAL_KIND), INTENT(IN) :: y
        TYPE(rpe_var) :: z
        z%sbits = MAX(significand_bits(x), significand_bits(y))
        z = x%val - y
    END FUNCTION sub_rpe_real

    ELEMENTAL FUNCTION sub_rpe_realalt (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        REAL(KIND=RPE_ALTERNATE_KIND), INTENT(IN) :: y
        TYPE(rpe_var) :: z
        z%sbits = MAX(significand_bits(x), significand_bits(y))
        z = x%val - y
    END FUNCTION sub_rpe_realalt

    ELEMENTAL FUNCTION sub_integer_rpe (x, y) RESULT (z)
        INTEGER, INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        TYPE(rpe_var) :: z
        z%sbits = MAX(significand_bits(x), significand_bits(y))
        z = x - y%val
    END FUNCTION sub_integer_rpe

    ELEMENTAL FUNCTION sub_long_rpe (x, y) RESULT (z)
        INTEGER(KIND=8), INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        TYPE(rpe_var) :: z
        z%sbits = MAX(significand_bits(x), significand_bits(y))
        z = x - y%val
    END FUNCTION sub_long_rpe

    ELEMENTAL FUNCTION sub_real_rpe (x, y) RESULT (z)
        REAL(KIND=RPE_REAL_KIND), INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        TYPE(rpe_var) :: z
        z%sbits = MAX(significand_bits(x), significand_bits(y))
        z = x - y%val
    END FUNCTION sub_real_rpe

    ELEMENTAL FUNCTION sub_realalt_rpe (x, y) RESULT (z)
        REAL(KIND=RPE_ALTERNATE_KIND), INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        TYPE(rpe_var) :: z
        z%sbits = MAX(significand_bits(x), significand_bits(y))
        z = x - y%val
    END FUNCTION sub_realalt_rpe

    !-------------------------------------------------------------------
    ! Overloaded definitions for (*):
    !

    ELEMENTAL FUNCTION mul_rpe_rpe (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        TYPE(rpe_var) :: z
        z%sbits = MAX(significand_bits(x), significand_bits(y))
        z = x%val * y%val
    END FUNCTION mul_rpe_rpe

    ELEMENTAL FUNCTION mul_rpe_integer (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        INTEGER, INTENT(IN) :: y
        TYPE(rpe_var) :: z
        z%sbits = MAX(significand_bits(x), significand_bits(y))
        z = x%val * y
    END FUNCTION mul_rpe_integer

    ELEMENTAL FUNCTION mul_rpe_long (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        INTEGER(KIND=8), INTENT(IN) :: y
        TYPE(rpe_var) :: z
        z%sbits = MAX(significand_bits(x), significand_bits(y))
        z = x%val * y
    END FUNCTION mul_rpe_long

    ELEMENTAL FUNCTION mul_rpe_real (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        REAL(KIND=RPE_REAL_KIND), INTENT(IN) :: y
        TYPE(rpe_var) :: z
        z%sbits = MAX(significand_bits(x), significand_bits(y))
        z = x%val * y
    END FUNCTION mul_rpe_real

    ELEMENTAL FUNCTION mul_rpe_realalt (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        REAL(KIND=RPE_ALTERNATE_KIND), INTENT(IN) :: y
        TYPE(rpe_var) :: z
        z%sbits = MAX(significand_bits(x), significand_bits(y))
        z = x%val * y
    END FUNCTION mul_rpe_realalt

    ELEMENTAL FUNCTION mul_integer_rpe (x, y) RESULT (z)
        INTEGER, INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        TYPE(rpe_var) :: z
        z%sbits = MAX(significand_bits(x), significand_bits(y))
        z = x * y%val
    END FUNCTION mul_integer_rpe

    ELEMENTAL FUNCTION mul_long_rpe (x, y) RESULT (z)
        INTEGER(KIND=8), INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        TYPE(rpe_var) :: z
        z%sbits = MAX(significand_bits(x), significand_bits(y))
        z = x * y%val
    END FUNCTION mul_long_rpe

    ELEMENTAL FUNCTION mul_real_rpe (x, y) RESULT (z)
        REAL(KIND=RPE_REAL_KIND), INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        TYPE(rpe_var) :: z
        z%sbits = MAX(significand_bits(x), significand_bits(y))
        z = x * y%val
    END FUNCTION mul_real_rpe

    ELEMENTAL FUNCTION mul_realalt_rpe (x, y) RESULT (z)
        REAL(KIND=RPE_ALTERNATE_KIND), INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        TYPE(rpe_var) :: z
        z%sbits = MAX(significand_bits(x), significand_bits(y))
        z = x * y%val
    END FUNCTION mul_realalt_rpe

    !-------------------------------------------------------------------
    ! Overloaded definitions for (/):
    !

    ELEMENTAL FUNCTION div_rpe_rpe (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        TYPE(rpe_var) :: z
        z%sbits = MAX(significand_bits(x), significand_bits(y))
        z = x%val / y%val
    END FUNCTION div_rpe_rpe

    ELEMENTAL FUNCTION div_rpe_integer (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        INTEGER, INTENT(IN) :: y
        TYPE(rpe_var) :: z
        z%sbits = MAX(significand_bits(x), significand_bits(y))
        z = x%val / y
    END FUNCTION div_rpe_integer

    ELEMENTAL FUNCTION div_rpe_long (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        INTEGER(KIND=8), INTENT(IN) :: y
        TYPE(rpe_var) :: z
        z%sbits = MAX(significand_bits(x), significand_bits(y))
        z = x%val / y
    END FUNCTION div_rpe_long

    ELEMENTAL FUNCTION div_rpe_real (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        REAL(KIND=RPE_REAL_KIND), INTENT(IN) :: y
        TYPE(rpe_var) :: z
        z%sbits = MAX(significand_bits(x), significand_bits(y))
        z = x%val / y
    END FUNCTION div_rpe_real

    ELEMENTAL FUNCTION div_rpe_realalt (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        REAL(KIND=RPE_ALTERNATE_KIND), INTENT(IN) :: y
        TYPE(rpe_var) :: z
        z%sbits = MAX(significand_bits(x), significand_bits(y))
        z = x%val / y
    END FUNCTION div_rpe_realalt

    ELEMENTAL FUNCTION div_integer_rpe (x, y) RESULT (z)
        INTEGER, INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        TYPE(rpe_var) :: z
        z%sbits = MAX(significand_bits(x), significand_bits(y))
        z = x / y%val
    END FUNCTION div_integer_rpe

    ELEMENTAL FUNCTION div_long_rpe (x, y) RESULT (z)
        INTEGER(KIND=8), INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        TYPE(rpe_var) :: z
        z%sbits = MAX(significand_bits(x), significand_bits(y))
        z = x / y%val
    END FUNCTION div_long_rpe

    ELEMENTAL FUNCTION div_real_rpe (x, y) RESULT (z)
        REAL(KIND=RPE_REAL_KIND), INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        TYPE(rpe_var) :: z
        z%sbits = MAX(significand_bits(x), significand_bits(y))
        z = x / y%val
    END FUNCTION div_real_rpe

    ELEMENTAL FUNCTION div_realalt_rpe (x, y) RESULT (z)
        REAL(KIND=RPE_ALTERNATE_KIND), INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        TYPE(rpe_var) :: z
        z%sbits = MAX(significand_bits(x), significand_bits(y))
        z = x / y%val
    END FUNCTION div_realalt_rpe

    !-------------------------------------------------------------------
    ! Overloaded definitions for (.GE.):
    !

    ELEMENTAL FUNCTION ge_rpe_rpe (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        LOGICAL :: z
        z = x%val .GE. y%val
    END FUNCTION ge_rpe_rpe

    ELEMENTAL FUNCTION ge_rpe_integer (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        INTEGER, INTENT(IN) :: y
        LOGICAL :: z
        z = x%val .GE. y
    END FUNCTION ge_rpe_integer

    ELEMENTAL FUNCTION ge_rpe_long (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        INTEGER(KIND=8), INTENT(IN) :: y
        LOGICAL :: z
        z = x%val .GE. y
    END FUNCTION ge_rpe_long

    ELEMENTAL FUNCTION ge_rpe_real (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        REAL(KIND=RPE_REAL_KIND), INTENT(IN) :: y
        LOGICAL :: z
        z = x%val .GE. y
    END FUNCTION ge_rpe_real

    ELEMENTAL FUNCTION ge_rpe_realalt (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        REAL(KIND=RPE_ALTERNATE_KIND), INTENT(IN) :: y
        LOGICAL :: z
        z = x%val .GE. y
    END FUNCTION ge_rpe_realalt

    ELEMENTAL FUNCTION ge_integer_rpe (x, y) RESULT (z)
        INTEGER, INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        LOGICAL :: z
        z = x .GE. y%val
    END FUNCTION ge_integer_rpe

    ELEMENTAL FUNCTION ge_long_rpe (x, y) RESULT (z)
        INTEGER(KIND=8), INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        LOGICAL :: z
        z = x .GE. y%val
    END FUNCTION ge_long_rpe

    ELEMENTAL FUNCTION ge_real_rpe (x, y) RESULT (z)
        REAL(KIND=RPE_REAL_KIND), INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        LOGICAL :: z
        z = x .GE. y%val
    END FUNCTION ge_real_rpe

    ELEMENTAL FUNCTION ge_realalt_rpe (x, y) RESULT (z)
        REAL(KIND=RPE_ALTERNATE_KIND), INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        LOGICAL :: z
        z = x .GE. y%val
    END FUNCTION ge_realalt_rpe

    !-------------------------------------------------------------------
    ! Overloaded definitions for (.LE.):
    !

    ELEMENTAL FUNCTION le_rpe_rpe (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        LOGICAL :: z
        z = x%val .LE. y%val
    END FUNCTION le_rpe_rpe

    ELEMENTAL FUNCTION le_rpe_integer (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        INTEGER, INTENT(IN) :: y
        LOGICAL :: z
        z = x%val .LE. y
    END FUNCTION le_rpe_integer

    ELEMENTAL FUNCTION le_rpe_long (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        INTEGER(KIND=8), INTENT(IN) :: y
        LOGICAL :: z
        z = x%val .LE. y
    END FUNCTION le_rpe_long

    ELEMENTAL FUNCTION le_rpe_real (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        REAL(KIND=RPE_REAL_KIND), INTENT(IN) :: y
        LOGICAL :: z
        z = x%val .LE. y
    END FUNCTION le_rpe_real

    ELEMENTAL FUNCTION le_rpe_realalt (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        REAL(KIND=RPE_ALTERNATE_KIND), INTENT(IN) :: y
        LOGICAL :: z
        z = x%val .LE. y
    END FUNCTION le_rpe_realalt

    ELEMENTAL FUNCTION le_integer_rpe (x, y) RESULT (z)
        INTEGER, INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        LOGICAL :: z
        z = x .LE. y%val
    END FUNCTION le_integer_rpe

    ELEMENTAL FUNCTION le_long_rpe (x, y) RESULT (z)
        INTEGER(KIND=8), INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        LOGICAL :: z
        z = x .LE. y%val
    END FUNCTION le_long_rpe

    ELEMENTAL FUNCTION le_real_rpe (x, y) RESULT (z)
        REAL(KIND=RPE_REAL_KIND), INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        LOGICAL :: z
        z = x .LE. y%val
    END FUNCTION le_real_rpe

    ELEMENTAL FUNCTION le_realalt_rpe (x, y) RESULT (z)
        REAL(KIND=RPE_ALTERNATE_KIND), INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        LOGICAL :: z
        z = x .LE. y%val
    END FUNCTION le_realalt_rpe

    !-------------------------------------------------------------------
    ! Overloaded definitions for (.GT.):
    !

    ELEMENTAL FUNCTION gt_rpe_rpe (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        LOGICAL :: z
        z = x%val .GT. y%val
    END FUNCTION gt_rpe_rpe

    ELEMENTAL FUNCTION gt_rpe_integer (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        INTEGER, INTENT(IN) :: y
        LOGICAL :: z
        z = x%val .GT. y
    END FUNCTION gt_rpe_integer

    ELEMENTAL FUNCTION gt_rpe_long (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        INTEGER(KIND=8), INTENT(IN) :: y
        LOGICAL :: z
        z = x%val .GT. y
    END FUNCTION gt_rpe_long

    ELEMENTAL FUNCTION gt_rpe_real (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        REAL(KIND=RPE_REAL_KIND), INTENT(IN) :: y
        LOGICAL :: z
        z = x%val .GT. y
    END FUNCTION gt_rpe_real

    ELEMENTAL FUNCTION gt_rpe_realalt (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        REAL(KIND=RPE_ALTERNATE_KIND), INTENT(IN) :: y
        LOGICAL :: z
        z = x%val .GT. y
    END FUNCTION gt_rpe_realalt

    ELEMENTAL FUNCTION gt_integer_rpe (x, y) RESULT (z)
        INTEGER, INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        LOGICAL :: z
        z = x .GT. y%val
    END FUNCTION gt_integer_rpe

    ELEMENTAL FUNCTION gt_long_rpe (x, y) RESULT (z)
        INTEGER(KIND=8), INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        LOGICAL :: z
        z = x .GT. y%val
    END FUNCTION gt_long_rpe

    ELEMENTAL FUNCTION gt_real_rpe (x, y) RESULT (z)
        REAL(KIND=RPE_REAL_KIND), INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        LOGICAL :: z
        z = x .GT. y%val
    END FUNCTION gt_real_rpe

    ELEMENTAL FUNCTION gt_realalt_rpe (x, y) RESULT (z)
        REAL(KIND=RPE_ALTERNATE_KIND), INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        LOGICAL :: z
        z = x .GT. y%val
    END FUNCTION gt_realalt_rpe

    !-------------------------------------------------------------------
    ! Overloaded definitions for (.LT.):
    !

    ELEMENTAL FUNCTION lt_rpe_rpe (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        LOGICAL :: z
        z = x%val .LT. y%val
    END FUNCTION lt_rpe_rpe

    ELEMENTAL FUNCTION lt_rpe_integer (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        INTEGER, INTENT(IN) :: y
        LOGICAL :: z
        z = x%val .LT. y
    END FUNCTION lt_rpe_integer

    ELEMENTAL FUNCTION lt_rpe_long (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        INTEGER(KIND=8), INTENT(IN) :: y
        LOGICAL :: z
        z = x%val .LT. y
    END FUNCTION lt_rpe_long

    ELEMENTAL FUNCTION lt_rpe_real (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        REAL(KIND=RPE_REAL_KIND), INTENT(IN) :: y
        LOGICAL :: z
        z = x%val .LT. y
    END FUNCTION lt_rpe_real

    ELEMENTAL FUNCTION lt_rpe_realalt (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        REAL(KIND=RPE_ALTERNATE_KIND), INTENT(IN) :: y
        LOGICAL :: z
        z = x%val .LT. y
    END FUNCTION lt_rpe_realalt

    ELEMENTAL FUNCTION lt_integer_rpe (x, y) RESULT (z)
        INTEGER, INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        LOGICAL :: z
        z = x .LT. y%val
    END FUNCTION lt_integer_rpe

    ELEMENTAL FUNCTION lt_long_rpe (x, y) RESULT (z)
        INTEGER(KIND=8), INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        LOGICAL :: z
        z = x .LT. y%val
    END FUNCTION lt_long_rpe

    ELEMENTAL FUNCTION lt_real_rpe (x, y) RESULT (z)
        REAL(KIND=RPE_REAL_KIND), INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        LOGICAL :: z
        z = x .LT. y%val
    END FUNCTION lt_real_rpe

    ELEMENTAL FUNCTION lt_realalt_rpe (x, y) RESULT (z)
        REAL(KIND=RPE_ALTERNATE_KIND), INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        LOGICAL :: z
        z = x .LT. y%val
    END FUNCTION lt_realalt_rpe

    !-------------------------------------------------------------------
    ! Overloaded definitions for (==):
    !

    ELEMENTAL FUNCTION eq_rpe_rpe (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        LOGICAL :: z
        z = x%val == y%val
    END FUNCTION eq_rpe_rpe

    ELEMENTAL FUNCTION eq_rpe_integer (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        INTEGER, INTENT(IN) :: y
        LOGICAL :: z
        z = x%val == y
    END FUNCTION eq_rpe_integer

    ELEMENTAL FUNCTION eq_rpe_long (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        INTEGER(KIND=8), INTENT(IN) :: y
        LOGICAL :: z
        z = x%val == y
    END FUNCTION eq_rpe_long

    ELEMENTAL FUNCTION eq_rpe_real (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        REAL(KIND=RPE_REAL_KIND), INTENT(IN) :: y
        LOGICAL :: z
        z = x%val == y
    END FUNCTION eq_rpe_real

    ELEMENTAL FUNCTION eq_rpe_realalt (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        REAL(KIND=RPE_ALTERNATE_KIND), INTENT(IN) :: y
        LOGICAL :: z
        z = x%val == y
    END FUNCTION eq_rpe_realalt

    ELEMENTAL FUNCTION eq_integer_rpe (x, y) RESULT (z)
        INTEGER, INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        LOGICAL :: z
        z = x == y%val
    END FUNCTION eq_integer_rpe

    ELEMENTAL FUNCTION eq_long_rpe (x, y) RESULT (z)
        INTEGER(KIND=8), INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        LOGICAL :: z
        z = x == y%val
    END FUNCTION eq_long_rpe

    ELEMENTAL FUNCTION eq_real_rpe (x, y) RESULT (z)
        REAL(KIND=RPE_REAL_KIND), INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        LOGICAL :: z
        z = x == y%val
    END FUNCTION eq_real_rpe

    ELEMENTAL FUNCTION eq_realalt_rpe (x, y) RESULT (z)
        REAL(KIND=RPE_ALTERNATE_KIND), INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        LOGICAL :: z
        z = x == y%val
    END FUNCTION eq_realalt_rpe

    !-------------------------------------------------------------------
    ! Overloaded definitions for (/=):
    !

    ELEMENTAL FUNCTION ne_rpe_rpe (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        LOGICAL :: z
        z = x%val /= y%val
    END FUNCTION ne_rpe_rpe

    ELEMENTAL FUNCTION ne_rpe_integer (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        INTEGER, INTENT(IN) :: y
        LOGICAL :: z
        z = x%val /= y
    END FUNCTION ne_rpe_integer

    ELEMENTAL FUNCTION ne_rpe_long (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        INTEGER(KIND=8), INTENT(IN) :: y
        LOGICAL :: z
        z = x%val /= y
    END FUNCTION ne_rpe_long

    ELEMENTAL FUNCTION ne_rpe_real (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        REAL(KIND=RPE_REAL_KIND), INTENT(IN) :: y
        LOGICAL :: z
        z = x%val /= y
    END FUNCTION ne_rpe_real

    ELEMENTAL FUNCTION ne_rpe_realalt (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        REAL(KIND=RPE_ALTERNATE_KIND), INTENT(IN) :: y
        LOGICAL :: z
        z = x%val /= y
    END FUNCTION ne_rpe_realalt

    ELEMENTAL FUNCTION ne_integer_rpe (x, y) RESULT (z)
        INTEGER, INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        LOGICAL :: z
        z = x /= y%val
    END FUNCTION ne_integer_rpe

    ELEMENTAL FUNCTION ne_long_rpe (x, y) RESULT (z)
        INTEGER(KIND=8), INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        LOGICAL :: z
        z = x /= y%val
    END FUNCTION ne_long_rpe

    ELEMENTAL FUNCTION ne_real_rpe (x, y) RESULT (z)
        REAL(KIND=RPE_REAL_KIND), INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        LOGICAL :: z
        z = x /= y%val
    END FUNCTION ne_real_rpe

    ELEMENTAL FUNCTION ne_realalt_rpe (x, y) RESULT (z)
        REAL(KIND=RPE_ALTERNATE_KIND), INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        LOGICAL :: z
        z = x /= y%val
    END FUNCTION ne_realalt_rpe

    !-------------------------------------------------------------------
    ! Overloaded definitions for (**):
    !

    ELEMENTAL FUNCTION pow_rpe_rpe (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        TYPE(rpe_var) :: z
        z%sbits = MAX(significand_bits(x), significand_bits(y))
        z = x%val ** y%val
    END FUNCTION pow_rpe_rpe

    ELEMENTAL FUNCTION pow_rpe_integer (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        INTEGER, INTENT(IN) :: y
        TYPE(rpe_var) :: z
        z%sbits = MAX(significand_bits(x), significand_bits(y))
        z = x%val ** y
    END FUNCTION pow_rpe_integer

    ELEMENTAL FUNCTION pow_rpe_long (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        INTEGER(KIND=8), INTENT(IN) :: y
        TYPE(rpe_var) :: z
        z%sbits = MAX(significand_bits(x), significand_bits(y))
        z = x%val ** y
    END FUNCTION pow_rpe_long

    ELEMENTAL FUNCTION pow_rpe_real (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        REAL(KIND=RPE_REAL_KIND), INTENT(IN) :: y
        TYPE(rpe_var) :: z
        z%sbits = MAX(significand_bits(x), significand_bits(y))
        z = x%val ** y
    END FUNCTION pow_rpe_real

    ELEMENTAL FUNCTION pow_rpe_realalt (x, y) RESULT (z)
        TYPE(rpe_var), INTENT(IN) :: x
        REAL(KIND=RPE_ALTERNATE_KIND), INTENT(IN) :: y
        TYPE(rpe_var) :: z
        z%sbits = MAX(significand_bits(x), significand_bits(y))
        z = x%val ** y
    END FUNCTION pow_rpe_realalt

    ELEMENTAL FUNCTION pow_integer_rpe (x, y) RESULT (z)
        INTEGER, INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        TYPE(rpe_var) :: z
        z%sbits = MAX(significand_bits(x), significand_bits(y))
        z = x ** y%val
    END FUNCTION pow_integer_rpe

    ELEMENTAL FUNCTION pow_long_rpe (x, y) RESULT (z)
        INTEGER(KIND=8), INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        TYPE(rpe_var) :: z
        z%sbits = MAX(significand_bits(x), significand_bits(y))
        z = x ** y%val
    END FUNCTION pow_long_rpe

    ELEMENTAL FUNCTION pow_real_rpe (x, y) RESULT (z)
        REAL(KIND=RPE_REAL_KIND), INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        TYPE(rpe_var) :: z
        z%sbits = MAX(significand_bits(x), significand_bits(y))
        z = x ** y%val
    END FUNCTION pow_real_rpe

    ELEMENTAL FUNCTION pow_realalt_rpe (x, y) RESULT (z)
        REAL(KIND=RPE_ALTERNATE_KIND), INTENT(IN) :: x
        TYPE(rpe_var), INTENT(IN) :: y
        TYPE(rpe_var) :: z
        z%sbits = MAX(significand_bits(x), significand_bits(y))
        z = x ** y%val
    END FUNCTION pow_realalt_rpe