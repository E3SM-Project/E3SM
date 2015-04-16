#include "private.h"
#include <stdio.h>

static float meanhashvalue (Hashentry *, int);

void GPTLprint_hashstats (FILE *fp, int nthreads, Hashentry **hashtable, int tablesize)
{
  int t;                    /* thread index */
  int i, ii;
  int totent;               /* per-thread collision count (diagnostic) */
  int nument;               /* per-index collision count (diagnostic) */
  /*
  ** Diagnostics for collisions and GPTL memory usage
  */
  int num_zero;             /* number of buckets with 0 collisions */
  int num_one;              /* number of buckets with 1 collision */
  int num_two;              /* number of buckets with 2 collisions */
  int num_more;             /* number of buckets with more than 2 collisions */
  int most;                 /* biggest collision count */
  bool first;

  for (t = 0; t < nthreads; t++) {
    first = true;
    totent   = 0;
    num_zero = 0;
    num_one  = 0;
    num_two  = 0;
    num_more = 0;
    most     = 0;

    for (i = 0; i < tablesize; i++) {
      nument = hashtable[t][i].nument;
      if (nument > 1) {
	totent += nument-1;
	if (first) {
	  first = false;
	  fprintf (fp, "\nthread %d had some hash collisions:\n", t);
	}
	fprintf (fp, "hashtable[%d][%d] had %d entries:", t, i, nument);
	for (ii = 0; ii < nument; ii++)
	  fprintf (fp, " %s", hashtable[t][i].entries[ii]->name);
	fprintf (fp, "\n");
      }
      switch (nument) {
      case 0:
	++num_zero;
	break;
      case 1:
	++num_one;
	break;
      case 2:
	++num_two;
	break;
      default:
	++num_more;
	break;
      }
      most = MAX (most, nument);
    }
    
    if (totent > 0) {
      fprintf (fp, "Total collisions thread %d = %d\n", t, totent);
      fprintf (fp, "Entry information:\n");
      fprintf (fp, "num_zero = %d num_one = %d num_two = %d num_more = %d\n",
	       num_zero, num_one, num_two, num_more);
      fprintf (fp, "Most = %d\n", most);
    }
  }
  fprintf (fp, "Size of hash table was %d\n", tablesize);
  fprintf (fp, "Mean hash index for thread 0 was %f\n", meanhashvalue (hashtable[0], tablesize));
}
  
static float meanhashvalue (Hashentry *hashtable, int tablesize)
{
  float sum = 0.;  /* used to calculate mean */
  int nument;
  int totent = 0;  /* number of entries */
  int i;
  
  for (i = 1; i < tablesize; ++i) {
    nument = hashtable[i].nument;
    if (nument > 0) {
      sum += (float) (nument * i);
      totent += hashtable[i].nument;
    }
  }
  if (totent == 0)
    return (float) 0.;
  else
    return sum / totent;
}
