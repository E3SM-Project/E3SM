    PUBLIC :: OPERATOR(+)
    INTERFACE OPERATOR(+)
        MODULE PROCEDURE add_rpe
        MODULE PROCEDURE add_rpe_rpe
        MODULE PROCEDURE add_rpe_integer
        MODULE PROCEDURE add_rpe_long
        MODULE PROCEDURE add_rpe_real
        MODULE PROCEDURE add_rpe_realalt
        MODULE PROCEDURE add_integer_rpe
        MODULE PROCEDURE add_long_rpe
        MODULE PROCEDURE add_real_rpe
        MODULE PROCEDURE add_realalt_rpe
    END INTERFACE OPERATOR(+)

    PUBLIC :: OPERATOR(-)
    INTERFACE OPERATOR(-)
        MODULE PROCEDURE sub_rpe
        MODULE PROCEDURE sub_rpe_rpe
        MODULE PROCEDURE sub_rpe_integer
        MODULE PROCEDURE sub_rpe_long
        MODULE PROCEDURE sub_rpe_real
        MODULE PROCEDURE sub_rpe_realalt
        MODULE PROCEDURE sub_integer_rpe
        MODULE PROCEDURE sub_long_rpe
        MODULE PROCEDURE sub_real_rpe
        MODULE PROCEDURE sub_realalt_rpe
    END INTERFACE OPERATOR(-)

    PUBLIC :: OPERATOR(*)
    INTERFACE OPERATOR(*)
        MODULE PROCEDURE mul_rpe_rpe
        MODULE PROCEDURE mul_rpe_integer
        MODULE PROCEDURE mul_rpe_long
        MODULE PROCEDURE mul_rpe_real
        MODULE PROCEDURE mul_rpe_realalt
        MODULE PROCEDURE mul_integer_rpe
        MODULE PROCEDURE mul_long_rpe
        MODULE PROCEDURE mul_real_rpe
        MODULE PROCEDURE mul_realalt_rpe
    END INTERFACE OPERATOR(*)

    PUBLIC :: OPERATOR(/)
    INTERFACE OPERATOR(/)
        MODULE PROCEDURE div_rpe_rpe
        MODULE PROCEDURE div_rpe_integer
        MODULE PROCEDURE div_rpe_long
        MODULE PROCEDURE div_rpe_real
        MODULE PROCEDURE div_rpe_realalt
        MODULE PROCEDURE div_integer_rpe
        MODULE PROCEDURE div_long_rpe
        MODULE PROCEDURE div_real_rpe
        MODULE PROCEDURE div_realalt_rpe
    END INTERFACE OPERATOR(/)

    PUBLIC :: OPERATOR(.GE.)
    INTERFACE OPERATOR(.GE.)
        MODULE PROCEDURE ge_rpe_rpe
        MODULE PROCEDURE ge_rpe_integer
        MODULE PROCEDURE ge_rpe_long
        MODULE PROCEDURE ge_rpe_real
        MODULE PROCEDURE ge_rpe_realalt
        MODULE PROCEDURE ge_integer_rpe
        MODULE PROCEDURE ge_long_rpe
        MODULE PROCEDURE ge_real_rpe
        MODULE PROCEDURE ge_realalt_rpe
    END INTERFACE OPERATOR(.GE.)

    PUBLIC :: OPERATOR(.LE.)
    INTERFACE OPERATOR(.LE.)
        MODULE PROCEDURE le_rpe_rpe
        MODULE PROCEDURE le_rpe_integer
        MODULE PROCEDURE le_rpe_long
        MODULE PROCEDURE le_rpe_real
        MODULE PROCEDURE le_rpe_realalt
        MODULE PROCEDURE le_integer_rpe
        MODULE PROCEDURE le_long_rpe
        MODULE PROCEDURE le_real_rpe
        MODULE PROCEDURE le_realalt_rpe
    END INTERFACE OPERATOR(.LE.)

    PUBLIC :: OPERATOR(.GT.)
    INTERFACE OPERATOR(.GT.)
        MODULE PROCEDURE gt_rpe_rpe
        MODULE PROCEDURE gt_rpe_integer
        MODULE PROCEDURE gt_rpe_long
        MODULE PROCEDURE gt_rpe_real
        MODULE PROCEDURE gt_rpe_realalt
        MODULE PROCEDURE gt_integer_rpe
        MODULE PROCEDURE gt_long_rpe
        MODULE PROCEDURE gt_real_rpe
        MODULE PROCEDURE gt_realalt_rpe
    END INTERFACE OPERATOR(.GT.)

    PUBLIC :: OPERATOR(.LT.)
    INTERFACE OPERATOR(.LT.)
        MODULE PROCEDURE lt_rpe_rpe
        MODULE PROCEDURE lt_rpe_integer
        MODULE PROCEDURE lt_rpe_long
        MODULE PROCEDURE lt_rpe_real
        MODULE PROCEDURE lt_rpe_realalt
        MODULE PROCEDURE lt_integer_rpe
        MODULE PROCEDURE lt_long_rpe
        MODULE PROCEDURE lt_real_rpe
        MODULE PROCEDURE lt_realalt_rpe
    END INTERFACE OPERATOR(.LT.)

    PUBLIC :: OPERATOR(==)
    INTERFACE OPERATOR(==)
        MODULE PROCEDURE eq_rpe_rpe
        MODULE PROCEDURE eq_rpe_integer
        MODULE PROCEDURE eq_rpe_long
        MODULE PROCEDURE eq_rpe_real
        MODULE PROCEDURE eq_rpe_realalt
        MODULE PROCEDURE eq_integer_rpe
        MODULE PROCEDURE eq_long_rpe
        MODULE PROCEDURE eq_real_rpe
        MODULE PROCEDURE eq_realalt_rpe
    END INTERFACE OPERATOR(==)

    PUBLIC :: OPERATOR(/=)
    INTERFACE OPERATOR(/=)
        MODULE PROCEDURE ne_rpe_rpe
        MODULE PROCEDURE ne_rpe_integer
        MODULE PROCEDURE ne_rpe_long
        MODULE PROCEDURE ne_rpe_real
        MODULE PROCEDURE ne_rpe_realalt
        MODULE PROCEDURE ne_integer_rpe
        MODULE PROCEDURE ne_long_rpe
        MODULE PROCEDURE ne_real_rpe
        MODULE PROCEDURE ne_realalt_rpe
    END INTERFACE OPERATOR(/=)

    PUBLIC :: OPERATOR(**)
    INTERFACE OPERATOR(**)
        MODULE PROCEDURE pow_rpe_rpe
        MODULE PROCEDURE pow_rpe_integer
        MODULE PROCEDURE pow_rpe_long
        MODULE PROCEDURE pow_rpe_real
        MODULE PROCEDURE pow_rpe_realalt
        MODULE PROCEDURE pow_integer_rpe
        MODULE PROCEDURE pow_long_rpe
        MODULE PROCEDURE pow_real_rpe
        MODULE PROCEDURE pow_realalt_rpe
    END INTERFACE OPERATOR(**)