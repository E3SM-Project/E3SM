/* ------------------------------------------------------------------ */
/*    function get_memory()                                           */
/*                                                                    */
/*                  +-----------------------------------+             */
/*                  | WAVEWATCH III           NOAA/NCEP |             */
/*                  |           H. L. Tolman            |             */
/*                  |                                 C |             */
/*                  | Last update :         28-Jan-2014 |             */
/*                  +-----------------------------------+             */
/*                                                                    */
/*    28-Jan-2014 : Origination.                     ( version 5.00 ) */
/*                                                                    */
/* 1. Purpose :                                                       */
/*                                                                    */
/*    Get memory use HWM in Mb for profiling purposes.                */
/*                                                                    */
/* 2. Method :                                                        */
/*                                                                    */
/*    C system calls, based on code mostly form Jim Abeles, provided  */
/*    by George Vandenberghe.                                         */
/*                                                                    */
/* 3. Parameters :                                                    */
/*                                                                    */
/* 4. Subroutines used :                                              */
/*                                                                    */
/* 5. Called by :                                                     */
/*                                                                    */
/* 6. Error messages :                                                */
/*                                                                    */
/* 7. Remarks :                                                       */
/*                                                                    */
/* 8. Structure :                                                     */
/*                                                                    */
/*    See source code.                                                */
/*                                                                    */
/* 9. Source code :                                                   */
/*                                                                    */
/* ------------------------------------------------------------------ */

#include <stdio.h>
#include <sys/resource.h>

double get_memory_(void)

{
  long Kbytes;
  double Mbytes;
  struct rusage RU;

  getrusage(RUSAGE_SELF, &RU);
 
  Kbytes = RU.ru_maxrss;

  Mbytes = ((double) Kbytes) / 1024.0;

  return Mbytes;
}

/* End of get_memory ------------------------------------------------ */
