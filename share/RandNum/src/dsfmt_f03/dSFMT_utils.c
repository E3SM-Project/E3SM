#include "dSFMT.h"
#include <stdlib.h>

/* copyright James Spencer 2012.
 * New BSD License, see License.txt for details.
 */

/* Utility (memory-access) functions to enable use of dSFMT from Fortran. */

void* malloc_dsfmt_t(void) {
    /* Allocate sufficient memory for a dSFMT state (ie a variable of type dsfmt_t). */
    return malloc(sizeof(dsfmt_t));
}

void free_dsfmt_t(dsfmt_t* ptr) {
    /* Free memory associated with a dSFMT state (ie a variable of type dsfmt_t). */
    free(ptr);
}
