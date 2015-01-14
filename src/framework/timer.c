#include <stdlib.h>
#include <stdio.h>

#ifdef UNDERSCORE
#define start_timer start_timer_
#define stop_timer stop_timer_
#else
#ifdef DOUBLEUNDERSCORE
#define start_timer start_timer__
#define stop_timer stop_timer__
#endif
#endif

#ifdef GETTIMEOFDAY
#include <sys/time.h>
#endif

#ifdef DARWIN
#include <mach/mach.h>
#include <mach/mach_time.h>
#include <unistd.h>
#endif

#ifdef LINUX
#include <time.h>
#endif

#ifdef GETTIMEOFDAY
struct timeval start_time[10];
struct timeval stop_time[10];
#endif

#ifdef DARWIN
uint64_t start_time[10];
uint64_t end_time[10];
#endif

#ifdef AIX
timebasestruct_t start_time[10];
timebasestruct_t end_time[10];
#endif

#ifdef LINUX
struct timespec start_time[10];
struct timespec end_time[10];
#endif

void start_timer_(int * n)
{
#ifdef GETTIMEOFDAY
   gettimeofday(&start_time[*n], NULL);
#endif

#ifdef DARWIN
   start_time[*n] = mach_absolute_time();
#endif

#ifdef AIX
   read_real_time(&start_time[*n], TIMEBASE_SZ);
#endif

#ifdef LINUX
   clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start_time[*n]);
#endif
}

void stop_timer_(int * n, int * secs, int * n_secs)
{
#ifdef GETTIMEOFDAY
   gettimeofday(&stop_time[*n], NULL);
  
   *secs   = (int)(stop_time[*n].tv_sec - start_time[*n].tv_sec);
   *n_secs = (int)(stop_time[*n].tv_usec - start_time[*n].tv_usec) * 1000;

   if (*n_secs < 0)  {
      *secs   -= 1;
      *n_secs += 1000000000;
   }
#endif

#ifdef DARWIN
   uint64_t elapsed, elapsedNano;
   static mach_timebase_info_data_t sTimebaseInfo;

   end_time[*n] = mach_absolute_time();

   elapsed = end_time[*n] - start_time[*n];

    if ( sTimebaseInfo.denom == 0 ) {
        (void) mach_timebase_info(&sTimebaseInfo);
    }

    // Do the maths. We hope that the multiplication doesn't 
    // overflow; the price you pay for working in fixed point.

    elapsedNano = elapsed * sTimebaseInfo.numer / sTimebaseInfo.denom;


   *secs   = (int)(elapsedNano / 1000000000);
   *n_secs = (int)(elapsedNano % 1000000000);
#endif

#ifdef AIX
   read_real_time(&end_time[*n], TIMEBASE_SZ);
   time_base_to_time(&start_time[*n], TIMEBASE_SZ);
   time_base_to_time(&end_time[*n], TIMEBASE_SZ);

   *secs = end_time[*n].tb_high - start_time[*n].tb_high;
   *n_secs = end_time[*n].tb_low - start_time[*n].tb_low;

   if (*n_secs < 0)  {
      *secs   -= 1;
      *n_secs += 1000000000;
   }
#endif

#ifdef LINUX
   clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end_time[*n]);

   *secs = (int)(end_time[*n].tv_sec - start_time[*n].tv_sec);
   *n_secs = (int)(end_time[*n].tv_nsec - start_time[*n].tv_nsec);

   if (*n_secs < 0)  {
      *secs   -= 1;
      *n_secs += 1000000000;
   }
#endif
}

