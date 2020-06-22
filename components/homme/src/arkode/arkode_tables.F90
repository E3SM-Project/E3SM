#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#define ARK324_ARK 1
#define ARK436_ARK 2
#define ARK453_ARK 3
#define ARS222_ARK 4
#define ARS232_ARK 5
#define ARS233_ARK 6
#define ARS343_ARK 7
#define ARS443_ARK 8
#define SSP3333B_ARK 9
#define SSP3333C_ARK 10
#define RK2_ARK 11
#define KGU35_ARK 12
#define IMKG232_ARK 13
#define IMKG242_ARK 14
#define IMKG243_ARK 15
#define IMKG252_ARK 16
#define IMKG253_ARK 17
#define IMKG254_ARK 18
#define IMKG342_ARK 19
#define IMKG343_ARK 20
#define IMKG353_ARK 21
#define IMKG354_ARK 22

module arkode_tables

  use kinds, only: real_kind

  implicit none

  private

  integer, parameter :: max_stage_num = 10

  ! data type for passing ARKode Butcher table names
  type :: table_list
    integer :: ARK324   = ARK324_ARK
    integer :: ARK436   = ARK436_ARK
    integer :: ARK453   = ARK453_ARK
    integer :: ARS222   = ARS222_ARK
    integer :: ARS232   = ARS232_ARK
    integer :: ARS233   = ARS233_ARK
    integer :: ARS343   = ARS343_ARK
    integer :: ARS443   = ARS443_ARK
    integer :: SSP3333B = SSP3333B_ARK
    integer :: SSP3333C = SSP3333C_ARK
    integer :: RK2      = RK2_ARK
    integer :: KGU35    = KGU35_ARK
    integer :: IMKG232  = IMKG232_ARK
    integer :: IMKG242  = IMKG242_ARK
    integer :: IMKG243  = IMKG243_ARK
    integer :: IMKG252  = IMKG252_ARK
    integer :: IMKG253  = IMKG253_ARK
    integer :: IMKG254  = IMKG254_ARK
    integer :: IMKG342  = IMKG342_ARK
    integer :: IMKG343  = IMKG343_ARK
    integer :: IMKG353  = IMKG353_ARK
    integer :: IMKG354  = IMKG353_ARK
  end type table_list

  ! data type for specifying Butcher table(s)
  type :: butcher_table_set
    ! *RK Method Information
    integer         :: imex ! 0=implicit, 1=explicit, 2=imex
    integer         :: s ! number of stages
    integer         :: q ! method order
    integer         :: p ! embedded method order
    ! Explicit Butcher Table
    real(real_kind) :: Ae(max_stage_num,max_stage_num)
    real(real_kind) :: be(max_stage_num)
    real(real_kind) :: be2(max_stage_num)
    real(real_kind) :: ce(max_stage_num)
    ! Implicit Butcher Table
    real(real_kind) :: Ai(max_stage_num,max_stage_num)
    real(real_kind) :: bi(max_stage_num)
    real(real_kind) :: bi2(max_stage_num)
    real(real_kind) :: ci(max_stage_num)
  end type butcher_table_set

  public :: table_list, butcher_table_set, set_Butcher_tables

  save

contains

  subroutine set_Butcher_tables(table_set, name)
    !-----------------------------------------------------------------
    ! Description: sets Butcher tables for ARKode solver
    !   Arguments:
    !     table_set - (butcher_table, in/output) object for butcher table(s)
    !          name - (integer, input) constant identifying table name
    !-----------------------------------------------------------------

    !======= Inclusions ===========
    use kinds,            only: real_kind
    use parallel_mod,     only: abortmp

    !======= Declarations =========
    implicit none

    ! calling variables
    type(butcher_table_set), target, intent(inout) :: table_set
    integer,                         intent(in)    :: name

    ! local variables
    type(butcher_table_set), pointer :: ts
    real(real_kind) :: beta, gamma, delta, b1, b2
    real(real_kind) :: a(max_stage_num), ahat(max_stage_num)
    real(real_kind) :: b(max_stage_num), dhat(max_stage_num)

    !======= Internals ============
    ts => table_set

    select case (name)

    case (RK2_ARK)
        ts%imex = 1 ! explicit
        ts%s = 2 ! 2 stage
        ts%q = 2 ! 2nd order
        ts%p = 0 ! no embedded order
        ts%be2 = 0.d0 ! no embedded explicit method
        ! Explicit Butcher Table (matrix)
        ts%Ae(1:2,1:2) = 0.d0
        ts%Ae(2,1) = 0.5d0
        ! Explicit Butcher Table (vectors)
        ts%ce(1:2) = (/ 0.d0, 0.5d0 /)
        ts%be(1:2) = (/ 0.d0, 1.d0 /)

      case (KGU35_ARK)
        ts%imex = 1 ! explicit
        ts%s = 5 ! 5 stage
        ts%q = 3 ! 3rd order
        ts%p = 0 ! no embedded order
        ts%be2 = 0.d0 ! no embedded explicit method
        ! Explicit Butcher Table (matrix)
        ts%Ae(1:5,1:5) = 0.d0
        ts%Ae(2,1) = 0.2d0
        ts%Ae(3,1:2) = (/  0.d0, 0.2d0 /)
        ts%Ae(4,1:3) = (/  0.d0,  0.d0, 1.d0/3.d0 /)
        ts%Ae(5,1:4) = (/  0.d0,  0.d0,      0.d0, 2.d0/3.d0 /)
        ! Explicit Butcher Table (vectors)
        ts%ce(1:5) = (/   0.d0, 0.2d0, 0.2d0, 1.d0/3.d0, 2.d0/3.d0 /)
        ts%be(1:5) = (/ 0.25d0,  0.d0,  0.d0,      0.d0,    0.75d0 /)

      case (ARS232_ARK)
        ts%imex = 2 ! imex
        ts%s = 3 ! 3 stage
        ts%q = 2 ! 2nd order
        ts%p = 0 ! no embedded order
        ts%be2 = 0.d0 ! no embedded explicit method
        delta = -2.d0*sqrt(2.d0)/3.d0
        gamma = 1.d0 - 1.d0/sqrt(2.d0)
        ! Implicit Butcher Table (matrix)
        ts%Ai(1:3,1:3) = 0.d0
        ts%Ai(2,1:2) = (/ 0.d0,      gamma /)
        ts%Ai(3,1:3) = (/ 0.d0, 1.d0-gamma, gamma /)
        ! Implicit Butcher Table (vectors)
        ts%ci(1:3) = (/ 0.d0, gamma, 1.d0 /)
        ts%bi(1:3) = (/ 0.d0, 1.d0-gamma, gamma /)
        ! Explicit Butcher Table (matrix)
        ts%Ae(1:3,1:3) = 0.d0
        ts%Ae(2,1) =  gamma
        ts%Ae(3,1:2) = (/ delta, 1.d0-delta /)
        ! Explicit Butcher Table (vectors)
        ts%ce(1:3) = ts%ci(1:3)
        ts%be(1:3) = ts%bi(1:3)

      case (ARK453_ARK)
        ts%imex = 2 ! imex
        ts%s = 5 ! 5 stage
        ts%q = 3 ! 3rd order
        ts%p = 0 ! no embedded order
        ts%be2 = 0.d0 ! no embedded explicit method
        ! Implicit Butcher Table (matrix)
        ts%Ai(1:5,1:5) = 0.d0
        ts%Ai(2,1:2) = (/ -0.22284985318525410d0, 0.32591194130117247d0 /)
        ts%Ai(3,1:3) = (/ -0.46801347074080545d0, 0.86349284225716961d0, &
                            0.32591194130117247d0 /)
        ts%Ai(4,1:4) = (/ -0.46509906651927421d0, 0.81063103116959553d0, &
                            0.61036726756832357d0, 0.32591194130117247d0 /)
        ts%Ai(5,1:5) = (/ 0.87795339639076675d0, -0.72692641526151547d0, &
                            0.75204137157372720d0, -0.22898029400415088d0, &
                            0.32591194130117247d0 /)
        ! Implicit Butcher Table (vectors)
        ts%ci(1:5) = (/ 0.d0, 0.1030620881159184d0, &
                          0.72139131281753662d0, 1.28181117351981733d0, &
                          1.d0 /)
        ts%bi(1:5) = (/ 0.87795339639076672d0, -0.72692641526151549d0, &
                          0.7520413715737272d0, -0.22898029400415090d0, &
                          0.32591194130117246d0 /)
        ! Explicit Butcher Table (matrix)
        ts%Ae(1:5,1:5) = 0.d0
        ts%Ae(2,1) = 0.10306208811591838d0
        ts%Ae(3,1:2) = (/ -0.94124866143519894d0, 1.6626399742527356d0 /)
        ts%Ae(4,1:3) = (/ -1.3670975201437765d0, 1.3815852911016873d0, &
                            1.2673234025619065d0 /)
        ts%Ae(5,1:4) = (/ -0.81287582068772448d0, 0.81223739060505738d0, &
                            0.90644429603699305d0, 0.094194134045674111d0 /)
        ! Explicit Butcher Table (vectors)
        ts%ce(1:5) = ts%ci(1:5)
        ts%be(1:5) = ts%bi(1:5)

      case (ARS233_ARK)
        ts%imex = 2 ! imex
        ts%s = 3 ! 3 stage
        ts%q = 3 ! 3rd order
        ts%p = 0 ! no embedded order
        ts%be2 = 0.d0 ! no embedded explicit method
        gamma = 1.d0/6.d0*(3.d0 + sqrt(3.d0))
        ! Implicit Butcher Table (matrix)
        ts%Ai(1:3,1:3) = 0.d0
        ts%Ai(2,1:2) = (/ 0.d0,          gamma /)
        ts%Ai(3,1:3) = (/ 0.d0, 1.d0-2.d0*gamma, gamma /)
        ! Implicit Butcher Table (vectors)
        ts%ci(1:3) = (/ 0.d0, gamma, 1.d0-gamma /)
        ts%bi(1:3) = (/ 0.d0, 0.5d0,      0.5d0 /)
        ! Explicit Butcher Table (matrix)
        ts%Ae(1:3,1:3) = 0.d0
        ts%Ae(2,1) = gamma
        ts%Ae(3,1:2) = (/ gamma-1.d0, 2.d0*(1.d0-gamma) /)
        ! Explicit Butcher Table (vectors)
        ts%ce(1:3) = ts%ci(1:3)
        ts%be(1:3) = ts%bi(1:3)

      case (ARS222_ARK)
        ts%imex = 2 ! imex
        ts%s = 3 ! 3 stage
        ts%q = 2 ! 2rd order
        ts%p = 0 ! no embedded order
        ts%be2 = 0.d0 ! no embedded explicit method
        gamma = 1.d0 - 1.d0/sqrt(2.d0)
        delta = 1.d0 - 1.d0/(2.d0*gamma)
        ! Implicit Butcher Table (matrix)
        ts%Ai(1:3,1:3) = 0.d0
        ts%Ai(2,1:2) = (/ 0.d0,      gamma /)
        ts%Ai(3,1:3) = (/ 0.d0, 1.d0-gamma, gamma /)
        ! Implicit Butcher Table (vectors)
        ts%ci(1:3) = (/ 0.d0,      gamma,  1.d0 /)
        ts%bi(1:3) = (/ 0.d0, 1.d0-gamma, gamma /)
        ! Explicit Butcher Table (matrix)
        ts%Ae(1:3,1:3) = 0.d0
        ts%Ae(2,1) = gamma
        ts%Ae(3,1:2) = (/ delta, 1.d0-delta /)
        ! Explicit Butcher Table (vectors)
        ts%ce(1:3) = ts%ci(1:3)
        ts%be(1:3) = ts%bi(1:3)

      case (ARS343_ARK)
        ts%imex = 2 ! imex
        ts%s = 4 ! 4 stage
        ts%q = 3 ! 3rd order
        ts%p = 0 ! no embedded order
        ts%be2 = 0.d0 ! no embedded explicit method
        gamma = 0.4358665215084590d0
        b1 = -3.d0*gamma*gamma/2.d0 + 4.d0*gamma - 1.d0/4.d0
        b2 = 3.d0*gamma*gamma/2.d0 - 5.d0*gamma + 5.d0/4.d0
        ! Implicit Butcher Table (matrix)
        ts%Ai(1:4,1:4) = 0.d0
        ts%Ai(2,1:2) = (/ 0.d0, gamma /)
        ts%Ai(3,1:3) = (/ 0.d0, (1.d0-gamma)/2.d0, gamma /)
        ts%Ai(4,1:4) = (/ 0.d0,                b1,    b2, gamma /)
        ! Implicit Butcher Table (vectors)
        ts%ci(1:4) = (/ 0.d0, gamma, (1.d0+gamma)/2.d0,  1.d0 /)
        ts%bi(1:4) = (/ 0.d0,    b1,                b2, gamma /)
        ! Explicit Butcher Table (matrix)
        ts%Ae(1:4,1:4) = 0.d0
        ts%Ae(2,1) = gamma
        ts%Ae(4,2) = 0.5529291480359398d0
        ts%Ae(4,3) = 0.5529291480359398d0
        ts%Ae(4,1) = 1.d0 - ts%Ae(4,2) - ts%Ae(4,3)
        ts%Ae(3,1) = ts%Ae(4,2)*(2.d0 - 9.d0*gamma + 3.d0*gamma*gamma)/2.d0 &
                    +ts%Ae(4,3)*(11.d0 - 42.d0*gamma + 15.d0*gamma*gamma)/4.d0 &
                    -7.d0/2.d0 + 13.d0*gamma - 9.d0*gamma*gamma/2.d0
        ts%Ae(3,2) = ts%Ae(4,2)*(-2.d0 + 9.d0*gamma - 3.d0*gamma*gamma)/2.d0 &
                    +ts%Ae(4,3)*(-11.d0 + 42.d0*gamma - 15.d0*gamma*gamma)/4.d0 &
                    +4.d0 - 25.d0*gamma/2.d0 + 9.d0*gamma*gamma/2.d0
        ! Explicit Butcher Table (vectors)
        ts%ce(1:4) = ts%ci(1:4)
        ts%be(1:4) = ts%bi(1:4)

      case (ARS443_ARK)
        ts%imex = 2 ! imex
        ts%s = 5 ! 5 stage
        ts%q = 3 ! 3rd order
        ts%p = 0 ! no embedded order
        ts%be2 = 0.d0 ! no embedded explicit method
        ! Implicit Butcher Table (matrix)
        ts%Ai(1:5,1:5) = 0.d0
        ts%Ai(2,1:2) = (/ 0.d0,  1.d0/2.d0 /)
        ts%Ai(3,1:3) = (/ 0.d0,  1.d0/6.d0,  1.d0/2.d0 /)
        ts%Ai(4,1:4) = (/ 0.d0, -1.d0/2.d0,  1.d0/2.d0, 1.d0/2.d0 /)
        ts%Ai(5,1:5) = (/ 0.d0,  3.d0/2.d0, -3.d0/2.d0, 1.d0/2.d0, 1.d0/2.d0 /)
        ! Implicit Butcher Table (vectors)
        ts%ci(1:5) = (/ 0.d0, 1.d0/2.d0,  2.d0/3.d0, 1.d0/2.d0,      1.d0 /)
        ts%bi(1:5) = (/ 0.d0, 3.d0/2.d0, -3.d0/2.d0, 1.d0/2.d0, 1.d0/2.d0 /)
        ! Explicit Butcher Table (matrix)
        ts%Ae(1:5,1:5) = 0.d0
        ts%Ae(2,1) = 1.d0/2.d0
        ts%Ae(3,1:2) = (/ 11.d0/18.d0, 1.d0/18.d0 /)
        ts%Ae(4,1:3) = (/   5.d0/6.d0, -5.d0/6.d0, 1.d0/2.d0 /)
        ts%Ae(5,1:4) = (/   1.d0/4.d0,  7.d0/4.d0, 3.d0/4.d0, -7.d0/4.d0 /)
        ! Explicit Butcher Table (vectors)
        ts%ce(1:5) = ts%ci(1:5)
        ts%be(1:5) = (/ 1.d0/4.d0, 7.d0/4.d0, 3.d0/4.d0, -7.d0/4.d0, 0.d0 /)

      case (ARK324_ARK)
        ts%imex = 2 ! imex
        ts%s = 4 ! 4 stage
        ts%q = 3 ! 3rd order
        ts%p = 2 ! 2nd order embedding
        ! Implicit Butcher Table (matrix)
        ts%Ai(1:4,1:4) = 0.d0
        ts%Ai(2,1:2) = (/ 1767732205903.d0/4055673282236.d0, &
                          1767732205903.d0/4055673282236.d0 /)
        ts%Ai(3,1:3) = (/ 2746238789719.d0/10658868560708.d0, &
                         -640167445237.d0/6845629431997.d0, &
                          1767732205903.d0/4055673282236.d0 /)
        ts%Ai(4,1:4) = (/ 1471266399579.d0/7840856788654.d0, &
                         -4482444167858.d0/7529755066697.d0, &
                          11266239266428.d0/11593286722821.d0, &
                          1767732205903.d0/4055673282236.d0 /)
        ! Implicit Butcher Table (vectors)
        ts%ci(1:4) = (/ 0.d0, 1767732205903.d0/2027836641118.d0, 3.d0/5.d0, 1.d0 /)
        ts%bi(1) = 1471266399579.d0/7840856788654.d0
        ts%bi(2) = -4482444167858.d0/7529755066697.d0
        ts%bi(3) =  11266239266428.d0/11593286722821.d0
        ts%bi(4) = 1767732205903.d0/4055673282236.d0
        ! Explicit Butcher Table (matrix)
        ts%Ae(1:4,1:4) = 0.d0
        ts%Ae(2,1) = 1767732205903.d0/2027836641118.d0
        ts%Ae(3,1:2) = (/ 5535828885825.d0/10492691773637.d0, &
                          788022342437.d0/10882634858940.d0 /)
        ts%Ae(4,1:3) = (/ 6485989280629.d0/16251701735622.d0, &
                         -4246266847089.d0/9704473918619.d0, &
                          10755448449292.d0/10357097424841.d0 /)
        ! Explicit Butcher Table (vectors)
        ts%ce(1:4) = ts%ci(1:4)
        ts%be(1:4) = ts%bi(1:4)
        ! Embedding
        ts%be2(1) = 2756255671327.d0/12835298489170.d0
        ts%be2(2) = -10771552573575.d0/22201958757719.d0
        ts%be2(3) = 9247589265047.d0/10645013368117.d0
        ts%be2(4) = 2193209047091.d0/5459859503100.d0

      case (ARK436_ARK)
        ts%imex = 2 ! imex
        ts%s = 6 ! 6 stage
        ts%q = 4 ! 4th order
        ts%p = 3 ! 3rd order embedding
        ! Implicit Butcher Table (matrix)
        ts%Ai(1:6,1:6) = 0.d0
        ts%Ai(2,1:2) = (/ 1.d0/4.d0, 1.d0/4.d0 /)
        ts%Ai(3,1:3) = (/ 8611.d0/62500.d0, -1743.d0/31250.d0, 1.d0/4.d0 /)
        ts%Ai(4,1:4) = (/ 5012029.d0/34652500.d0, -654441.d0/2922500.d0, &
                          174375.d0/388108.d0, 1.d0/4.d0 /)
        ts%Ai(5,1:5) = (/ 15267082809.d0/155376265600.d0, &
                          -71443401.d0/120774400.d0, 730878875.d0/902184768.d0, &
                          2285395.d0/8070912.d0, 1.d0/4.d0 /)
        ts%Ai(6,1:6) = (/ 82889.d0/524892.d0, 0.d0, 15625.d0/83664.d0, &
                          69875.d0/102672.d0, -2260.d0/8211.d0, 1.d0/4.d0 /)
        ! Implicit Butcher Table (vectors)
        ts%ci(1:6) = (/ 0.d0, 1.d0/2.d0, 83.d0/250.d0, 31.d0/50.d0, &
                      17.d0/20.d0, 1.d0 /)
        ts%bi(1:6) = (/ 82889.d0/524892.d0, 0.d0, 15625.d0/83664.d0, &
                      69875.d0/102672.d0, -2260.d0/8211.d0, 1.d0/4.d0 /)
        ! Explicit Butcher Table (matrix)
        ts%Ae(1:6,1:6) = 0.d0
        ts%Ae(2,1) = 1.d0/2.d0
        ts%Ae(3,1:2) = (/ 13861.d0/62500.d0, 6889.d0/62500.d0 /)
        ts%Ae(4,1:3) = (/ -116923316275.d0/2393684061468.d0, &
                          -2731218467317.d0/15368042101831.d0, &
                          9408046702089.d0/11113171139209.d0 /)
        ts%Ae(5,1:4) = (/ -451086348788.d0/2902428689909.d0, &
                          -2682348792572.d0/7519795681897.d0, &
                          12662868775082.d0/11960479115383.d0, &
                          3355817975965.d0/11060851509271.d0 /)
        ts%Ae(6,1:5) = (/ 647845179188.d0/3216320057751.d0, &
                          73281519250.d0/8382639484533.d0, &
                          552539513391.d0/3454668386233.d0, &
                          3354512671639.d0/8306763924573.d0, 4040.d0/17871.d0 /)
        ! Explicit Butcher Table (vectors)
        ts%ce(1:6) = ts%ci(1:6)
        ts%be(1:6) = ts%bi(1:6)
        ! Embedding
        ts%be2(1:6) = (/ 4586570599.d0/29645900160.d0, 0.d0, 178811875.d0/945068544.d0, &
                        814220225.d0/1159782912.d0, -3700637.d0/11593932.d0, &
                        61727.d0/225920.d0 /)

      case (SSP3333B_ARK)
        ts%imex = 2 ! imex
        ts%s = 3 ! 3 stage
        ts%q = 3 ! 3rd order
        ts%p = 0 ! no embedded order
        ts%be2 = 0.d0 ! no embedded explicit method
        ! Implicit Butcher Table (matrix)
        ts%Ai(1:3,1:3) = 0.d0
        ts%Ai(2,1:2) = (/      0.d0,       1.d0 /)
        ts%Ai(3,1:3) = (/ 1.d0/6.d0, -1.d0/3.d0, 2.d0/3.d0 /)
        ! Implicit Butcher Table (vectors)
        ts%ci(1:3) = (/ 0.d0, 1.d0, 0.5d0 /)
        ts%bi(1:3) = (/ 1.d0/6.d0, 1.d0/6.d0, 2.d0/3.d0 /)
        ! Explicit Butcher Table (matrix)
        ts%Ae(1:3,1:3) = 0.d0
        ts%Ae(2,1) = 1.d0
        ts%Ae(3,1:2) = (/ 0.25d0, 0.25d0 /)
        ! Explicit Butcher Table (vectors)
        ts%ce(1:3) = ts%ci(1:3)
        ts%be(1:3) = ts%bi(1:3)

      case (SSP3333C_ARK)
        ts%imex = 2 ! imex
        ts%s = 3 ! 3 stage
        ts%q = 3 ! 3rd order
        ts%p = 0 ! no embedded order
        ts%be2 = 0.d0 ! no embedded explicit method
        beta = sqrt(3.d0)/5.d0 + 0.5d0
        gamma = -1.d0/8.d0*(sqrt(3.d0) + 1.d0)
        ! Implicit Butcher Table (matrix)
        ts%Ai(1:3,1:3) = 0.d0
        ts%Ai(2,1:2) = (/ 4.d0*gamma + 2.d0*beta, 1.d0 - 4.d0*gamma - 2.d0*beta /)
        ts%Ai(3,1:3) = (/   0.5d0 - beta - gamma,                         gamma, beta /)
        ! Implicit Butcher Table (vectors)
        ts%ci(1:3) = (/ 0.d0, 1.d0, 0.5d0 /)
        ts%bi(1:3) = (/ 1.d0/6.d0, 1.d0/6.d0, 2.d0/3.d0 /)
        ! Explicit Butcher Table (matrix)
        ts%Ae(1:3,1:3) = 0.d0
        ts%Ae(2,1) = 1.d0
        ts%Ae(3,1:2) = (/ 0.25d0, 0.25d0 /)
        ! Explicit Butcher Table (vectors)
        ts%ce(1:3) = ts%ci(1:3)
        ts%be(1:3) = ts%bi(1:3)

      case (IMKG232_ARK)
        ts%imex = 2 ! imex
        ts%s = 4 ! 4 stage
        ts%q = 2 ! 2nd order
        ts%p = 0 ! no embedded order
        ts%be2 = 0.d0 ! no embedded explicit method
        ! IMEX-KG vectors
        delta = 0.5d0*(2.d0+sqrt(2.d0))
        a(1:3) = (/ 0.5d0, 0.5d0, 1.d0 /)
        ahat(1:3) = (/ 0.d0, -0.5d0*(sqrt(2.d0)+1.d0), 1.d0 /)
        dhat(1:2) = (/ delta, delta /)
        b(1:2) = 0.d0
        ! set IMEX-KG Butcher table
        call set_IMKG_Butcher_tables(table_set, ts%s, a, ahat, dhat, b)

      case (IMKG242_ARK)
        ts%imex = 2 ! imex
        ts%s = 5 ! 5 stage
        ts%q = 2 ! 2nd order
        ts%p = 0 ! no embedded order
        ts%be2 = 0.d0 ! no embedded explicit method
        ! IMEX-KG vectors
        delta = 0.5d0*(2.d0+sqrt(2.d0))
        a(1:4) = (/ 0.25d0, 1.d0/3.d0, 0.5d0, 1.d0 /)
        ahat(1:4) = (/ 0.d0, 0.d0, -0.5d0*(sqrt(2.d0)+1.d0), 1.d0 /)
        dhat(1:3) = (/ 0.d0, delta, delta /)
        b(1:3) = 0.d0
        ! set IMEX-KG Butcher table
        call set_IMKG_Butcher_tables(table_set, ts%s, a, ahat, dhat, b)

      case (IMKG243_ARK)
        ts%imex = 2 ! imex
        ts%s = 5 ! 5 stage
        ts%q = 2 ! 2nd order
        ts%p = 0 ! no embedded order
        ts%be2 = 0.d0 ! no embedded explicit method
        ! IMEX-KG vectors
        delta = 0.5d0+sqrt(3.d0)/6.d0
        a(1:4) = (/ 0.25d0, 1.d0/3.d0, 0.5d0, 1.d0 /)
        ahat(1:4) = (/ 0.d0, 1.d0/6.d0, -sqrt(3.d0)/6.d0, 1.d0 /)
        dhat(1:3) = (/ delta, delta, delta /)
        b(1:3) = 0.d0
        ! set IMEX-KG Butcher table
        call set_IMKG_Butcher_tables(table_set, ts%s, a, ahat, dhat, b)

      case (IMKG252_ARK)
        ts%imex = 2 ! imex
        ts%s = 6 ! 6 stage
        ts%q = 2 ! 2nd order
        ts%p = 0 ! no embedded order
        ts%be2 = 0.d0 ! no embedded explicit method
        ! IMEX-KG vectors
        delta = 0.5d0*(2.d0+sqrt(2.0))
        a(1:5) = (/ 0.25d0, 1.d0/6.d0, 3.d0/8.d0, 0.5d0, 1.d0 /)
        ahat(1:5) = (/ 0.d0, 0.d0, 0.d0, -0.5d0*(sqrt(2.d0)+1.d0), 1.d0 /)
        dhat(1:4) = (/ 0.d0, 0.d0, delta, delta /)
        b(1:4) = 0.d0
        ! set IMEX-KG Butcher table
        call set_IMKG_Butcher_tables(table_set, ts%s, a, ahat, dhat, b)

      case (IMKG253_ARK)
        ts%imex = 2 ! imex
        ts%s = 6 ! 6 stage
        ts%q = 2 ! 2nd order
        ts%p = 0 ! no embedded order
        ts%be2 = 0.d0 ! no embedded explicit method
        ! IMEX-KG vectors
        delta = 0.5d0+sqrt(3.d0)/6.d0
        gamma = 0.25d0*sqrt(3.d0)*(1.d0+sqrt(3.d0)/3.d0)*((sqrt(3.d0)/3.d0-1.d0)**2-2.d0)
        a(1:5) = (/ 0.25d0, 1.d0/6.d0, 3.d0/8.d0, 0.5d0, 1.d0 /)
        ahat(1:5) = (/ 0.d0, 0.d0, gamma, -sqrt(3.d0)/6.d0, 1.d0 /)
        dhat(1:4) = (/ 0.d0, delta, delta, delta /)
        b(1:4) = 0.d0
        ! set IMEX-KG Butcher table
        call set_IMKG_Butcher_tables(table_set, ts%s, a, ahat, dhat, b)

      case (IMKG254_ARK)
        ts%imex = 2 ! imex
        ts%s = 6 ! 6 stage
        ts%q = 2 ! 2nd order
        ts%p = 0 ! no embedded order
        ts%be2 = 0.d0 ! no embedded explicit method
        ! IMEX-KG vectors
        a(1:5) = (/ 0.25d0, 1.d0/6.d0, 3.d0/8.d0, 0.5d0, 1.d0 /)
        ahat(1:5) = (/ 0.d0, -1.d0/20.d0, 1.25d0, -0.5d0, 1.d0 /)
        dhat(1:4) = (/ -0.5d0, 1.d0, 1.d0, 1.d0 /)
        b(1:4) = 0.d0
        ! set IMEX-KG Butcher table
        call set_IMKG_Butcher_tables(table_set, ts%s, a, ahat, dhat, b)

      case (IMKG342_ARK)
        ts%imex = 2 ! imex
        ts%s = 5 ! 5 stage
        ts%q = 3 ! 3rd order
        ts%p = 0 ! no embedded order
        ts%be2 = 0.d0 ! no embedded explicit method
        ! IMEX-KG vectors
        delta = 0.5d0*(1.d0+sqrt(3.d0)/3.d0)
        a(1:4) = (/ 0.25d0, 2.d0/3.d0, 1.d0/3.d0, 0.75d0 /)
        ahat(1:4) = (/ 0.d0, (1.d0-sqrt(3.d0))/6.d0, -(1.d0+sqrt(3.d0))/6.d0, 0.75d0 /)
        dhat(1:3) = (/ 0.d0, delta, delta /)
        b(1:3) = (/ 0.d0, 1.d0/3.d0, 0.25d0 /)
        ! set IMEX-KG Butcher table
        call set_IMKG_Butcher_tables(table_set, ts%s, a, ahat, dhat, b)


      case (IMKG343_ARK)
        ts%imex = 2 ! imex
        ts%s = 5 ! 5 stage
        ts%q = 3 ! 3rd order
        ts%p = 0 ! no embedded order
        ts%be2 = 0.d0 ! no embedded explicit method
        ! IMEX-KG vectors
        a(1:4) = (/ 0.25d0, 2.d0/3.d0, 1.d0/3.d0, 0.75d0 /)
        ahat(1:4) = (/ 0.d0, -1.d0/3.d0, -2.d0/3.d0, 0.75d0 /)
        dhat(1:3) = (/ -1.d0/3.d0, 1.d0, 1.d0 /)
        b(1:3) = (/ 0.d0, 1.d0/3.d0, 0.25d0 /)
        ! set IMEX-KG Butcher table
        call set_IMKG_Butcher_tables(table_set, ts%s, a, ahat, dhat, b)

      case (IMKG353_ARK)
        ts%imex = 2 ! imex
        ts%s = 6 ! 6 stage
        ts%q = 3 ! 3rd order
        ts%p = 0 ! no embedded order
        ts%be2 = 0.d0 ! no embedded explicit method
        ! IMEX-KG vectors
        a(1:5) = (/ -0.017391304347826087d0, -.92d0, 5.d0/3.d0, 1.d0/3.d0, 3.d0/4.d0 /)
        ahat(1:2) = (/ 0.3075640504095504d0, -1.2990164859879263d0 /)
        ahat(3:5) = (/ 1.2516666666666665d0, -0.8166666666666668d0, 3d0/4d0/)
        dhat(1:4) = (/ -0.2981612530370581d0, .415d0, .415d0, 1.15d0 /)
        b(1:4) = (/ 1.d0, -1.d0, 1.d0/3.d0, 0.25d0 /)
        ! set IMEX-KG Butcher table
        call set_IMKG_Butcher_tables(table_set, ts%s, a, ahat, dhat, b)


      case default
        call abortmp('Unknown ARKode Butcher table name')
    end select

  end subroutine set_Butcher_tables
  !=================================================================

  subroutine set_IMKG_Butcher_tables(table_set, s, a, ahat, dhat, b)
    !-----------------------------------------------------------------
    ! Description: sets Butcher tables in IMEX-KG format:
    !   Ae:   0 |    0             Ai:     0 |       0
    !        c1 |   a1    0            chat1 |   ahat1 dhat1
    !        c2 |   b1   a2  0         chat2 |    b1 ahat2 dhat2
    !        c3 |   b2    0 a3 0       chat3 |    b2     0 ahat3 dhat3
    !           |                            |
    !        cq | bq-1 0 ... 0 aq 0    chatq |  bq-1     0    ...    0 ahatq 0
    !        ----------------------    -----------------------------------------
    !           | bq-1 0 ... 0 aq 0          |  bq-1     0    ...    0 ahatq 0
    !
    !        ci = sum_j Ae(i,j)        chati = sum_j Ai(i,j)
    !
    !   Arguments:
    !    arkode_parameters - (parameter, in/output) object for arkode parameters
    !                    s - (integer, input) number of stages
    !                    a - (real(max_stage_num), input) alpha vector
    !                 ahat - (real(max_stage_num), input) alpha_hat vector
    !                 dhat - (real(max_stage_num), in/output) delta_hat vector
    !                    b - (real(max_stage_num), input) beta vector
    !-----------------------------------------------------------------

    !======= Inclusions ===========
    use kinds,            only: real_kind
    use parallel_mod,     only: abortmp

    !======= Declarations =========
    implicit none

    ! calling variables
    type(butcher_table_set), target, intent(inout) :: table_set
    integer,                         intent(in)    :: s
    real(real_kind),                 intent(in)    :: a(max_stage_num)
    real(real_kind),                 intent(in)    :: ahat(max_stage_num)
    real(real_kind),                 intent(inout) :: dhat(max_stage_num)
    real(real_kind),                 intent(in)    :: b(max_stage_num)

    ! local variables
    type(butcher_table_set), pointer :: ts

    !======= Internals ============
    ts => table_set

    ! for easier implementation, extend dhat vector by setting s-1 index to zero
    dhat(s-1) = 0.d0

    ! set matrices and c vectors according to IMEX-KG format
    ts%Ai(1:s,1:s) = 0.d0
    ts%Ae(1:s,1:s) = 0.d0
    ts%ci(1) = 0.d0
    ts%ce(1) = 0.d0
    if (s > 1) then
      ts%Ai(2,1:2) = (/ ahat(1), dhat(1) /)
      ts%Ae(2,1) = a(1)
      ts%ci(2) = ahat(1)+dhat(1)
      ts%ce(2) = a(1)
    end if
    if (s > 2) then
      ts%Ai(3,1:3) = (/ b(1), ahat(2), dhat(2) /)
      ts%Ae(3,1:2) = (/ b(1), a(2) /)
      ts%ci(3) = b(1)+ahat(2)+dhat(2)
      ts%ce(3) = b(1)+a(2)
    end if
    if (s > 3) then
      ts%Ai(4,1:4) = (/ b(2), 0.d0, ahat(3), dhat(3) /)
      ts%Ae(4,1:3) = (/ b(2), 0.d0, a(3) /)
      ts%ci(4) = b(2)+ahat(3)+dhat(3)
      ts%ce(4) = b(2)+a(3)
    end if
    if (s > 4) then
      ts%Ai(5,1:5) = (/ b(3), 0.d0, 0.d0, ahat(4), dhat(4) /)
      ts%Ae(5,1:4) = (/ b(3), 0.d0, 0.d0, a(4) /)
      ts%ci(5) = b(3)+ahat(4)+dhat(4)
      ts%ce(5) = b(3)+a(4)
    end if
    if (s > 5) then
      ts%Ai(6,1:6) = (/ b(4), 0.d0, 0.d0, 0.d0, ahat(5), dhat(5) /)
      ts%Ae(6,1:5) = (/ b(4), 0.d0, 0.d0, 0.d0, a(5) /)
      ts%ci(6) = b(4)+ahat(5)+dhat(5)
      ts%ce(6) = b(4)+a(5)
    end if
    if (s > 6) then
      call abortmp('ARKode IMEX-KG only implemented for 6 stages or less')
    end if

    ! set b vectors
    ts%bi(1:s) = ts%Ai(s,1:s)
    ts%be(1:s) = ts%Ae(s,1:s)

  end subroutine set_IMKG_Butcher_tables

end module arkode_tables
