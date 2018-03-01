/*
 * src/x86.c
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2001
 *
 * Contains functions to set and restore the round-to-double flag in the
 * control word of a x86 FPU.
 */

#define _NO_CHANGE      0
#define _UPPER_CASE     1
#define _ADD_UNDERSCORE 2
#define _ADD_TWO_UNDERSCORES 3

#ifdef FORTRANUNDERSCORE
#define NAMING _ADD_UNDERSCORE
#endif

#ifdef FORTRANDOUBLEUNDERSCORE
#define NAMING _ADD_TWO_UNDERSCORES
#endif

#ifdef FORTRANCAPS
#define NAMING _UPPER_CASE
#endif

#ifndef NAMING
#define NAMING _NO_CHANGE
#endif

#if (NAMING == _ADD_UNDERSCORE)
#define x86_fix_start x86_fix_start_
#define x86_fix_end   x86_fix_end_
#endif

#if (NAMING == _ADD_TWO_UNDERSCORES)
#define x86_fix_start x86_fix_start__
#define x86_fix_end   x86_fix_end__
#endif

#if (NAMING == _UPPER_CASE)
#define x86_fix_start X86_FIX_START
#define x86_fix_end   X86_FIX_END
#endif

#ifdef x86
#ifndef _FPU_GETCW
#define _FPU_GETCW(x) asm volatile ("fnstcw %0":"=m" (x));
#endif

#ifndef _FPU_SETCW
#define _FPU_SETCW(x) asm volatile ("fldcw %0": :"m" (x));
#endif

#ifndef _FPU_EXTENDED
#define _FPU_EXTENDED 0x0300
#endif

#ifndef _FPU_DOUBLE
#define _FPU_DOUBLE 0x0200
#endif
#endif  /* x86 */

void x86_fix_start(unsigned short *old_cw) {
#ifdef x86
  unsigned short new_cw;

  _FPU_GETCW(*old_cw);
  new_cw = (*old_cw & ~_FPU_EXTENDED) | _FPU_DOUBLE;
  _FPU_SETCW(new_cw);
#endif
}

void x86_fix_end(unsigned short *old_cw) {
#ifdef x86
  _FPU_SETCW(*old_cw);
#endif
}

