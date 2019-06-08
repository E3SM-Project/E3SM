    PUBLIC :: epsilon
    INTERFACE epsilon
        MODULE PROCEDURE epsilon_rpe
    END INTERFACE epsilon

    PUBLIC :: tiny
    INTERFACE tiny
        MODULE PROCEDURE tiny_rpe
    END INTERFACE tiny

    PUBLIC :: kind
    INTERFACE kind
        MODULE PROCEDURE kind_rpe
    END INTERFACE kind

    PUBLIC :: abs
    INTERFACE abs
        MODULE PROCEDURE abs_rpe
    END INTERFACE abs

    PUBLIC :: cos
    INTERFACE cos
        MODULE PROCEDURE cos_rpe
    END INTERFACE cos

    PUBLIC :: sin
    INTERFACE sin
        MODULE PROCEDURE sin_rpe
    END INTERFACE sin

    PUBLIC :: tan
    INTERFACE tan
        MODULE PROCEDURE tan_rpe
    END INTERFACE tan

    PUBLIC :: acos
    INTERFACE acos
        MODULE PROCEDURE acos_rpe
    END INTERFACE acos

    PUBLIC :: asin
    INTERFACE asin
        MODULE PROCEDURE asin_rpe
    END INTERFACE asin

    PUBLIC :: atan
    INTERFACE atan
        MODULE PROCEDURE atan_rpe
    END INTERFACE atan

    PUBLIC :: cosh
    INTERFACE cosh
        MODULE PROCEDURE cosh_rpe
    END INTERFACE cosh

    PUBLIC :: sinh
    INTERFACE sinh
        MODULE PROCEDURE sinh_rpe
    END INTERFACE sinh

    PUBLIC :: tanh
    INTERFACE tanh
        MODULE PROCEDURE tanh_rpe
    END INTERFACE tanh

    PUBLIC :: exp
    INTERFACE exp
        MODULE PROCEDURE exp_rpe
    END INTERFACE exp

    PUBLIC :: log
    INTERFACE log
        MODULE PROCEDURE log_rpe
    END INTERFACE log

    PUBLIC :: log10
    INTERFACE log10
        MODULE PROCEDURE log10_rpe
    END INTERFACE log10

    PUBLIC :: sqrt
    INTERFACE sqrt
        MODULE PROCEDURE sqrt_rpe
    END INTERFACE sqrt

    PUBLIC :: spacing
    INTERFACE spacing
        MODULE PROCEDURE spacing_rpe
    END INTERFACE spacing

    PUBLIC :: floor
    INTERFACE floor
        MODULE PROCEDURE floor_rpe
    END INTERFACE floor

    PUBLIC :: int
    INTERFACE int
        MODULE PROCEDURE int_rpe
    END INTERFACE int

    PUBLIC :: nint
    INTERFACE nint
        MODULE PROCEDURE nint_rpe
    END INTERFACE nint

    PUBLIC :: atan2
    INTERFACE atan2
        MODULE PROCEDURE atan2_rpe_rpe
        MODULE PROCEDURE atan2_rpe_real
        MODULE PROCEDURE atan2_real_rpe
    END INTERFACE atan2

    PUBLIC :: dim
    INTERFACE dim
        MODULE PROCEDURE dim_rpe_rpe
        MODULE PROCEDURE dim_rpe_real
        MODULE PROCEDURE dim_real_rpe
    END INTERFACE dim

    PUBLIC :: mod
    INTERFACE mod
        MODULE PROCEDURE mod_rpe_rpe
        MODULE PROCEDURE mod_rpe_real
        MODULE PROCEDURE mod_real_rpe
    END INTERFACE mod

    PUBLIC :: nearest
    INTERFACE nearest
        MODULE PROCEDURE nearest_rpe_rpe
        MODULE PROCEDURE nearest_rpe_real
        MODULE PROCEDURE nearest_real_rpe
    END INTERFACE nearest

    PUBLIC :: sign
    INTERFACE sign
        MODULE PROCEDURE sign_rpe_rpe
        MODULE PROCEDURE sign_rpe_real
        MODULE PROCEDURE sign_real_rpe
    END INTERFACE sign

    PUBLIC :: min
    INTERFACE min
        MODULE PROCEDURE min_rpe_rpe
        MODULE PROCEDURE min_rpe_real
        MODULE PROCEDURE min_real_rpe
        MODULE PROCEDURE min_ma_rpe_rpe_rpe
        MODULE PROCEDURE min_ma_rpe_rpe_rpe_rpe
        MODULE PROCEDURE min_ma_rpe_rpe_rpe_rpe_rpe
        MODULE PROCEDURE min_ma_rpe_rpe_rpe_rpe_rpe_rpe
        MODULE PROCEDURE min_ma_rpe_rpe_rpe_rpe_rpe_rpe_rpe
        MODULE PROCEDURE min_ma_rpe_rpe_rpe_rpe_rpe_rpe_rpe_rpe
        MODULE PROCEDURE min_ma_rpe_rpe_rpe_rpe_rpe_rpe_rpe_rpe_rpe
    END INTERFACE min

    PUBLIC :: max
    INTERFACE max
        MODULE PROCEDURE max_rpe_rpe
        MODULE PROCEDURE max_rpe_real
        MODULE PROCEDURE max_real_rpe
        MODULE PROCEDURE max_ma_rpe_rpe_rpe
        MODULE PROCEDURE max_ma_rpe_rpe_rpe_rpe
        MODULE PROCEDURE max_ma_rpe_rpe_rpe_rpe_rpe
        MODULE PROCEDURE max_ma_rpe_rpe_rpe_rpe_rpe_rpe
        MODULE PROCEDURE max_ma_rpe_rpe_rpe_rpe_rpe_rpe_rpe
        MODULE PROCEDURE max_ma_rpe_rpe_rpe_rpe_rpe_rpe_rpe_rpe
        MODULE PROCEDURE max_ma_rpe_rpe_rpe_rpe_rpe_rpe_rpe_rpe_rpe
    END INTERFACE max

    PUBLIC :: minval
    INTERFACE minval
        MODULE PROCEDURE minval_rpe_1d
        MODULE PROCEDURE minval_rpe_2d
        MODULE PROCEDURE minval_rpe_3d
        MODULE PROCEDURE minval_rpe_4d
        MODULE PROCEDURE minval_rpe_5d
    END INTERFACE minval

    PUBLIC :: maxval
    INTERFACE maxval
        MODULE PROCEDURE maxval_rpe_1d
        MODULE PROCEDURE maxval_rpe_2d
        MODULE PROCEDURE maxval_rpe_3d
        MODULE PROCEDURE maxval_rpe_4d
        MODULE PROCEDURE maxval_rpe_5d
    END INTERFACE maxval

    PUBLIC :: sum
    INTERFACE sum
        MODULE PROCEDURE sum_rpe_1d
        MODULE PROCEDURE sum_rpe_2d
        MODULE PROCEDURE sum_rpe_3d
        MODULE PROCEDURE sum_rpe_4d
        MODULE PROCEDURE sum_rpe_5d
    END INTERFACE sum