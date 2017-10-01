module dynamic_vector_int_ptr

use ptr_wrapper, only: int_ptr

use pfunit_mod, only: throw

use shr_log_mod, only: OOBMsg => shr_log_OOBMsg

implicit none
private

#define VECTOR_NAME int_ptr_vector
#define TYPE_NAME type(int_ptr)
#define THROW(string) call throw(string)

public :: VECTOR_NAME

#include "dynamic_vector_typedef.inc"

contains

#include "dynamic_vector_procdef.inc"

end module dynamic_vector_int_ptr
