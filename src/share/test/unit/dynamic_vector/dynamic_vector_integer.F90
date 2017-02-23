module dynamic_vector_integer

use pfunit_mod, only: throw

use shr_log_mod, only: OOBMsg => shr_log_OOBMsg

implicit none
private

#define VECTOR_NAME integer_vector
#define TYPE_NAME integer
#define THROW(string) call throw(string)

public :: VECTOR_NAME

#include "dynamic_vector_typedef.inc"

contains

#include "dynamic_vector_procdef.inc"

end module dynamic_vector_integer
