module dynamic_vector_character16

use pfunit_mod, only: throw

use shr_log_mod, only: OOBMsg => shr_log_OOBMsg

implicit none
private

#define VECTOR_NAME character16_vector
#define TYPE_NAME character(len=16)
#define THROW(string) call throw(string)

public :: VECTOR_NAME

#include "dynamic_vector_typedef.inc"

contains

#include "dynamic_vector_procdef.inc"

end module dynamic_vector_character16
