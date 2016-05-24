/*----------------------------------------------------------*/
/* Simple elapsed-time timers for hand instrumenting codes. */
/*----------------------------------------------------------*/
/* link with libtimers.a                                    */
/*----------------------------------------------------------*/
/* Fortran:                                                 */
/*    call tbeg('label')  start timing                      */
/*    call tend('label')  stop  timing                      */
/*    call tprt()  print timer values and labels            */
/*----------------------------------------------------------*/
/* C / C++:                                                 */
/*    Tbeg("label");  start timing                          */
/*    Tend("label");  stop  timing                          */
/*    Tprt();  print timer values and labels                */
/*----------------------------------------------------------*/
/* author: Bob Walkup <walkup@us.ibm.com>                   */
/* modified by Theron Voran <voran@colorado.edu>            */
/*----------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/* #include <mpi.h> */

/*---------------------------------------*/
/* routine to read the timebase register */
/*---------------------------------------*/
#if defined(_AIX)
#include <sys/systemcfg.h>
void timebase(long long *);
#elif defined(BGL)
long long rts_get_timebase(void);
#else
#include <sys/time.h>
double timer_clock(void);
#endif


/*---------------------*/
/* Function prototypes */
/*---------------------*/
void Tbeg(char *);
void Tend(char *);
void Tprt(void);

void tbeg(char *, int);
void tend(char *, int);
void tprt(void);


static int index_from_label(char *);


/*---------------------------*/
/* variables with file-scope */
/*---------------------------*/
static int initialized = 0;
static int code_block  = 0;
static double tconv;

#define MAX_CODE_BLOCKS 500

static long long timer_in[MAX_CODE_BLOCKS];
static long long timer_sum[MAX_CODE_BLOCKS];
static char code_block_label[MAX_CODE_BLOCKS][80];
static int block_starts[MAX_CODE_BLOCKS];
static int block_stops[MAX_CODE_BLOCKS];



/*=======================================================*/
/* Fortran interface for tbeg: terminate the string      */
/* note: the length argument is hidden in Fortran source */
/*=======================================================*/
void tbeg(char * f_label, int length)
{
   int i, j, rc;
   double xint, xfrac;
   char this_label[80];

   strncpy(this_label, f_label, length);
   this_label[length] = '\0';

   /*---------------------------------------------------------------*/
   /* if the timers are not initialized, then do the initialization */
   /*---------------------------------------------------------------*/
   if (!initialized)
   {
       initialized = 1;

#if defined(_AIX)
       xint  = (double) _system_configuration.Xint;
       xfrac = (double) _system_configuration.Xfrac;
       tconv = 1.0e-9*xint/xfrac;
#elif defined(BGL)
       tconv = 1.0/700.0e6; 
#else
       tconv = 1.0/1.0e6;
#endif

       /*--------------------------------------*/
       /* set the initial timer values to zero */
       /*--------------------------------------*/
       for (j=0; j<MAX_CODE_BLOCKS; j++) timer_sum[j] = 0.0;

       /*------------------------------------------*/
       /* zero-out the code block starts and stops */
       /*------------------------------------------*/
       for (j=0; j<MAX_CODE_BLOCKS; j++)
       {
           block_starts[j] = 0;
           block_stops[j]  = 0;
       }

   }

   j = index_from_label(this_label);

   block_starts[j] += 1;

#if defined(_AIX)
   timebase(&timer_in[j]);
#elif defined(BGL)
   timer_in[j] = rts_get_timebase();
#else
   timer_in[j] = timer_clock();
#endif

   return;
}


/*=======================================================*/
/* Fortran interface for tend: terminate the string      */
/* note: the length argument is hidden in Fortran source */
/*=======================================================*/
void tend(char * f_label, int length)
{
   int i, j, rc;
   long long tnow;
   char this_label[80];

   strncpy(this_label, f_label, length);
   this_label[length] = '\0';

#if defined(_AIX)
   timebase(&tnow);
#elif defined(BGL)
   tnow = rts_get_timebase();
#else
   tnow = timer_clock();
#endif

   if (code_block >= MAX_CODE_BLOCKS) return;

   j = index_from_label(this_label);

   block_stops[j] += 1;

   timer_sum[j] += tnow - timer_in[j];

   return;
}

/*================================*/
/* Fortran interface for tprt     */
/*================================*/
void tprt(void)
{
   Tprt();
}


/*=====================================================*/
/* Initialize and start timing.                        */
/*=====================================================*/
void Tbeg(char * this_label)
{
   int i, j, rc;
   double xint, xfrac;

   /*---------------------------------------------------------------*/
   /* if the timers are not initialized, then do the initialization */
   /*---------------------------------------------------------------*/
   if (!initialized)
   {
       initialized = 1;

#if defined(_AIX)
       xint  = (double) _system_configuration.Xint;
       xfrac = (double) _system_configuration.Xfrac;
       tconv = 1.0e-9*xint/xfrac;
#elif defined(BGL)
       tconv = 1.0/700.0e6;
#else
       tconv = 1.0/1.0e6;
#endif

       /*--------------------------------------*/
       /* set the initial timer values to zero */
       /*--------------------------------------*/
       for (j=0; j<MAX_CODE_BLOCKS; j++)
              timer_sum[j] = 0.0;

       /*-------------------------------------------*/
       /* keep track of code block starts and stops */
       /*-------------------------------------------*/
       for (j=0; j<MAX_CODE_BLOCKS; j++)
       {
           block_starts[j] = 0;
           block_stops[j]  = 0;
       }
   }

   j = index_from_label(this_label);

   block_starts[j] += 1;

#if defined(_AIX)
   timebase(&timer_in[j]);
#elif defined(BGL)
   timer_in[j] = rts_get_timebase();
#else
   timer_in[j] = timer_clock();
#endif

   return;
}



/*================================================*/
/* stop timing, save the sum, and continue timing */
/*================================================*/
void Tend(char * this_label)
{
   int i, j, rc;
   long long tnow;

#if defined(_AIX)
   timebase(&tnow);
#elif defined(BGL)
   tnow = rts_get_timebase();
#else
   tnow = timer_clock();
#endif

   if (code_block >= MAX_CODE_BLOCKS) return;

   j = index_from_label(this_label);

   block_stops[j] += 1;

   timer_sum[j] += tnow - timer_in[j];

   return;
}



/*====================================*/
/* print the timer values with labels */
/*====================================*/
void Tprt(void)
{
   int rc, i, j, nblocks, taskid;
   double elapsed_seconds;
   char filename[80];
   char label[8][80];
   FILE * fp;

/*    MPI_Comm_rank(MPI_COMM_WORLD, &taskid); */
   taskid = 0;

   sprintf(filename, "timing_summary.%d", taskid);

   fp = fopen(filename, "w");
   if (fp == NULL)
   {
      fprintf(stderr, "from Tprt: could not open %s\n", filename);
      fprintf(stderr, "timing summary not printed\n");
      return;
   }

   fprintf(fp, "\n");
   fprintf(fp, "----------------------------------------------------------------\n");
   fprintf(fp, "Timing  summary:         #calls        time(sec)\n");
   fprintf(fp, "----------------------------------------------------------------\n");
   if (code_block >= MAX_CODE_BLOCKS) nblocks = MAX_CODE_BLOCKS;
   else                               nblocks = code_block;

   for (j=0; j<nblocks; j++)
   { 
       if (block_starts[j] == block_stops[j])
       {
           elapsed_seconds = tconv * ((double) timer_sum[j]);
           fprintf(fp, "%-20s  %9d  %14.3f\n", 
                  code_block_label[j], block_starts[j], elapsed_seconds);
       }
       else
       {
           fprintf(fp, "mismatch in starts/stops for code block '%s'\n", code_block_label[j]);
           fprintf(fp, "  starts = %d\n", block_starts[j]);
           fprintf(fp, "  stops  = %d\n", block_stops[j]);
       }
   }
   fprintf(fp, "\n");

   fclose(fp);

   return;
}


/*===========================================*/
/* Find the code-block number from the label.*/
/*===========================================*/
int index_from_label(char * this_label)
{
   int i, match;
   char * ptr;

   if (code_block < MAX_CODE_BLOCKS)
   {
      match = 0;
      for (i=code_block-1; i>=0; i--)
      {
         if (0 == strcmp(code_block_label[i], this_label))
         {
            match = 1;
            break;
         }
      }
    
      if (match == 0)
      {
         i = code_block;
         ptr = strcpy(code_block_label[i], this_label);
         if (ptr == NULL) code_block_label[i][0] = '\0';
         code_block ++;
      }
   }

   return i;

}


#if !defined(_AIX) && !defined(BGL)
double timer_clock()
{
  double t;
  struct timeval buffer;
  struct timezone dummy;
  gettimeofday (&buffer, &dummy);
  t = (double)(buffer.tv_sec*1000000 + buffer.tv_usec);
  return (t);
}
#endif

