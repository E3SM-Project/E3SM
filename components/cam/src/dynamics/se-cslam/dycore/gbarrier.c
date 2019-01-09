//: multi-level barrier code; predefined to a max of 64 threads, as below

// We need to define the Log2 of the maximum number of threads:
#define LOG2MAX 6
#define NTHREADS 64

#include <sched.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <omp.h>
#include <stdbool.h>
#include <stdint.h>

// utility functions:
int ipow2 (int val) {
  int result = 1;
  while (val > 0) {
    result *= 2;
    --val;
  }
  return result;
}

// Define the data associated with a global barrier:
typedef struct gbt {
  volatile bool LocalFlags    [2][LOG2MAX];
  volatile bool *PartnerFlags [2][LOG2MAX];
  bool sense;
  int  parity;
  int  id;
} GBarrier_Type;

// Define a singular type for the global barrier:
typedef struct gb {
  GBarrier_Type threadData[NTHREADS];
  int numThreads;
  int log2Threads;
} GBarrier;

void initializeThread(GBarrier_Type *threadData, int thread, int numThreads) {
  // Local loop variables: (p)arity, (r)ound and (x) [temporary]
  int p, r;
  unsigned int x;

  // local log2 threads:
  int log2Threads = ceil(log2(numThreads));

  threadData[thread].id = thread;
  threadData[thread].sense = true;
  threadData[thread].parity = 0;

  for (p = 0; p < 2; p++) {
    for (r = 0; r < log2Threads; r++) {
      x = (threadData[thread].id + ipow2(r)) % numThreads;
      threadData[thread].LocalFlags[p][r] = 0;
      threadData[thread].PartnerFlags[p][r] = &threadData[x].LocalFlags[p][r];
    }
  }
}

void gbarrier_synchronize(GBarrier* b, int thread)
{
  // Local:
  int i;

  // Get the pointer to our thread's data:
  GBarrier_Type *my = &b->threadData[thread];

  // Loop through the log2 rounds:
  for (i = 0; i < b->log2Threads; i++) {
    *my->PartnerFlags[my->parity][i] = my->sense;

    while (my->LocalFlags[my->parity][i] != my->sense) { sched_yield(); }
  }

  // Reverse the sense for reuse on parity=1
  if (my->parity == 1) { my->sense = !my->sense; }

  // Swap our parity between 0 & 1:
  my->parity = 1 - my->parity;
}

void gbarrier_initialize(GBarrier **ptb, int numThreads) {
  // Local variables:
  int t;

  GBarrier *b;
  (*ptb) = malloc(sizeof(GBarrier));
  b = (*ptb);

  b->numThreads = numThreads;
  b->log2Threads = ceil(log2(b->numThreads));

  for (t = 0; t < b->numThreads; t++) {
    initializeThread(b->threadData, t, b->numThreads);
  }
}

void gbarrier_print(GBarrier *b) {
  printf("GBarrier Info: %d threads \n", b->numThreads);
}

void gbarrier_free(GBarrier **ptb) {
  GBarrier *b = (*ptb);
  free(b);
}

