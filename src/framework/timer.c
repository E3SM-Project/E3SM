#include <stdlib.h>
#include <stdio.h>

#include <sys/time.h>

struct timeval start_time[10];
struct timeval stop_time[10];

void start_timer(int n)
{
   gettimeofday(&start_time[n], NULL);
}

void stop_timer(int n, int * secs, int * u_secs)
{
   gettimeofday(&stop_time[n], NULL);
  
   *secs   = (int)(stop_time[n].tv_sec - start_time[n].tv_sec);
   *u_secs = (int)(stop_time[n].tv_usec - start_time[n].tv_usec);

   if (*u_secs < 0)  {
      *secs   -= 1;
      *u_secs += 1000000;
   }
}

