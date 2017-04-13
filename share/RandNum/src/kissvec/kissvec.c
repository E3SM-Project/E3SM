// KISS random generator implemented in C
// Basic algorithm from George Marsaglia circa 1998
// Public domain Fortran implementation from http://www.fortran.com/
// downloaded by pjr on 03/16/04
// converted to vector form, functions inlined by pjr,mvr on 05/10/2004
// Translated back into C in 2015

#include <stddef.h>
#include <stdint.h>

#define shiftl_xor(x, n) (x ^= (x << n))
#define shiftr_xor(x, n) (x ^= (x >> n))

// The  KISS (Keep It Simple Stupid) random number generator. Combines:
// (1) The congruential generator x(n)=69069*x(n-1)+1327217885, period 2^32.
// (2) A 3-shift shift-register generator, period 2^32-1,
// (3) Two 16-bit multiply-with-carry generators, period 597273182964842497>2^59
//  Overall period>2^123;
//
// Note use of the C99 restrict keyword to enable optimization.
void kiss_rng(uint32_t seed1[restrict], uint32_t seed2[restrict],
              uint32_t seed3[restrict], uint32_t seed4[restrict],
              double ran_arr[restrict], size_t length) {
  size_t i;

  for (i = 0; i < length; ++i) {
    seed1[i] = 69069U * seed1[i] + 1327217885U;
    shiftl_xor(seed2[i], 13);
    shiftr_xor(seed2[i], 17);
    shiftl_xor(seed2[i], 5);
    seed3[i] = 18000U * (seed3[i] & 65535U) + (seed3[i] >> 16);
    seed4[i] = 30903U * (seed4[i] & 65535U) + (seed4[i] >> 16);
    ran_arr[i] = ((int32_t) (seed1[i] + seed2[i] + (seed3[i] << 16) + seed4[i])) * 2.328306E-10 + 0.5;
  }

}
