module dynamic_vector_r8

use shr_kind_mod, only: r8 => shr_kind_r8
use shr_infnan_mod, only: assignment(=), nan => shr_infnan_nan

use pfunit_mod, only: throw

use shr_log_mod, only: OOBMsg => shr_log_OOBMsg

implicit none
private

#define VECTOR_NAME r8_vector
#define TYPE_NAME real(r8)
#define THROW(string) call throw(string)

public :: VECTOR_NAME

#include "dynamic_vector_typedef.inc"

contains

#include "dynamic_vector_procdef.inc"

end module dynamic_vector_r8
