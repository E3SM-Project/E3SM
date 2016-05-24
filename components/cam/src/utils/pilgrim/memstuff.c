#if defined (USE_MLP)
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <errno.h>

#define CACHE_LINE_SIZE 128
#define SECTION_ROUND (1024*1024)
#define BARRIER_AREA_ROUND (16*1024)
#define round_up(n, round)  (((n + (round - 1)) / round) * round)

#if defined ( PGI )  || defined ( LAHEY )
#define ROUND (1024*1024)

volatile long long *__prior_counter_ptr;
volatile long long *__current_counter_ptr; 
#endif


/**********************************************************************/
/**********************************************************************/
/* Barrier                                                            */
/**********************************************************************/
#if defined ( PGI ) || defined ( LAHEY )


char *_lock_filename = "./_lock_file";
void
_lock(void)    /*This part will be implemented by semaphore later, BW */
{
  unsigned long long num_attempts = 0;
  int fd;

  while ((fd = creat(_lock_filename, 0)) == -1) {
    if (errno != EACCES) {
      perror("Failed to creat _lock_file");
      exit(1);
    }
/*    if ((++num_attempts % 100000) == 0) {       */
    if ((++num_attempts % 100000) == 0) {
      fprintf(stderr, "%lld lock attemps have failed on %s\n"
                "(Perhaps an old copy of the file was not deleted?)\n",
                num_attempts, _lock_filename);
    }
  }
  close(fd);
}

void
_unlock(void) /*This part will be implemented by semaphore later, BW */
{
  if (unlink(_lock_filename) == -1) {
    perror("unlink of _lock_file");
    exit(1);
  }
}                         

void
mlp_barrier__(int *n_ptr)    
{
  long long target_value = *n_ptr + *__prior_counter_ptr;
  long long new_value;

/*  printf("getpid = %d \n", getpid());   */
/*  printf("target_value = %d for %d \n", target_value, getpid()); */
  _lock();
  {
    new_value = *__current_counter_ptr + 1;
    *__current_counter_ptr = new_value;
  }
  _unlock();

  if (new_value == target_value) {
    *__prior_counter_ptr = target_value;
  }
  while (*__prior_counter_ptr != target_value) {
    /* spin */ ;
  }
}                                 

#else


/* Note: FAN_IN is arbitrary - it does not need to be a power of two.
** By using a power of two, the compiler can optimize the arithmetic.
** FAN_IN of 2 would be a binary tree, FAN_IN of 4 a quad tree, etc.
*/
#define FAN_IN 4
#define MAX_PROCESSES 1024


typedef struct {
    volatile unsigned long long current_value;
    unsigned long long pad1[CACHE_LINE_SIZE/sizeof(unsigned long long) - 1];
    volatile unsigned long long previous_value;
    unsigned long long pad2[CACHE_LINE_SIZE/sizeof(unsigned long long) - 1];
} __mlp_join_line_type;

__mlp_join_line_type *__ptr_mlp_join_area;





/*
** Barrier wait for n processes with id's  0 .. n-1
**
** Each block of FAN_IN processes increments a different counter.  When
** all the procs in that block have reported in, the group counter one
** level up is incremented.  When all the groups on that level have reported
** in, the next higher level counter is incremented.  When everyone has
** reported in, notification is propagated back down and we are done.
**
** Note: we must acquire the target_value *before* we __add_and_fetch
** the counter to avoid a race condition.
**
** By using 64bit monotonic counters, we do not need to reset the counts.
** This saves memory transactions.  The code as written does NOT assume
** FAN_IN is a power of two; we rely on the compiler to simplify the
** arithmetic where possible and appropriate.
** (we use unsigned values to make this more clear to the compiler).
*/



/* Recurse up the barrier tree */
/* Note that the id's are zero-based */
void
__mlp_sync_step(
    __mlp_join_line_type *base,
    unsigned int my_id,
    unsigned int last_id)
{
    unsigned int last_group = last_id / FAN_IN;
    unsigned int my_group = my_id / FAN_IN; 
    unsigned int my_group_size =
        (my_group < last_group) ?  FAN_IN : (last_id % FAN_IN) + 1;

    __mlp_join_line_type *ptr_my_group = base + my_group;

    unsigned long long target_value =
        ptr_my_group->previous_value + my_group_size;

    unsigned long long value =
        __add_and_fetch(&(ptr_my_group->current_value), 1ULL);

    if (value < target_value) {
        /* I am not the last process of my group that needs to get here */
        /* Spin until previous_value is set by the last member of the group */
        while (ptr_my_group->previous_value < target_value) /* spin */ ;
    }
    else {
        /* I was the last process in my group to reach this sync point. */

        /* If we are not already at the top level, then sync one level up */
        if (last_group != 0) {
          __mlp_sync_step(base + last_group + 1, my_group, last_group);
        }

        /* Propagate notification downwards to other processes in my group */
        ptr_my_group->previous_value = target_value;
    }
}


void
mlp_barrier_(unsigned int *ptr_my_id, unsigned int *ptr_num_procs)
{
    __mlp_sync_step(__ptr_mlp_join_area, *ptr_my_id, *ptr_num_procs - 1);
}

#endif

/**********************************************************************/
/**********************************************************************/
/* End Barrier                                                        */
/**********************************************************************/







/**********************************************************************/
/**********************************************************************/
/* Allocate and Initialize the shared memory                          */
/**********************************************************************/
void
mlp_getmem_(long long *n_ptr, long long size[], long long pointer[])
{
  int fd, i, n = *n_ptr;
  unsigned long long total_size = 0;
  char buf[100];
  char *mmap_addr;
#if defined ( IRIX64 )
  unsigned long barrier_area_size = 0;
  unsigned long num_groups_this_level = 0;
#endif

  for (i=0; i<n; i++) {
    total_size += size[i];
  }
        
  /* Add space for the barrier counters */

#if defined ( PGI ) || defined ( LAHEY )
  total_size += 2*sizeof(long long);
  total_size = ((total_size + (ROUND - 1))/ROUND) * ROUND;     
#else
  num_groups_this_level = round_up(MAX_PROCESSES, FAN_IN) / FAN_IN;
  for (;;) {
    barrier_area_size += sizeof(__mlp_join_line_type)*num_groups_this_level;
    if (num_groups_this_level == 1) break;
    num_groups_this_level = round_up(num_groups_this_level, FAN_IN) / FAN_IN;
  }

  barrier_area_size = round_up(barrier_area_size, BARRIER_AREA_ROUND);
  total_size += barrier_area_size;
        
  /* round up the whole thing */
  total_size = round_up(total_size, SECTION_ROUND);

#endif

  /* Generate and open a temp file */
  sprintf(buf, "./tempfile.%d", getpid());
  unlink(buf);
  fd = open(buf, O_RDWR|O_CREAT, 0600);
  if (fd < 0) {
    perror("open(2) of the shm file");
    exit(1);
  }
  /* Unlink the temp file so it will be automatically removed */
  unlink(buf);

  /* Grow the file */
  if (ftruncate(fd, (off_t)total_size) < 0) {
    perror("ftruncate(2) of the shm file");
    exit(1); 
  }
  
  /* Map the file "shared" */
  mmap_addr = (char *)mmap(0, (size_t) total_size,
                        PROT_READ|PROT_WRITE, MAP_SHARED, fd, 0);
  if (mmap_addr == (char *)MAP_FAILED) {
    perror("mmap(2) of the shm file");
    exit(1);
  }

#if defined ( IRIX64 )
  /* Initialize the barrier area */
  __ptr_mlp_join_area = (__mlp_join_line_type *) mmap_addr;
  for (i = barrier_area_size/sizeof(__mlp_join_line_type); i >=0; i--) {
    __ptr_mlp_join_area[i].current_value = 0;
    __ptr_mlp_join_area[i].previous_value = 0;
  }
  mmap_addr += barrier_area_size;
#endif

  /* Parcel out the space */
  for (i=0; i<n; i++) {
    pointer[i] = (long long) mmap_addr;
    mmap_addr += size[i];
  }

#if defined ( PGI ) || defined ( LAHEY )
  __prior_counter_ptr = (volatile long long *) mmap_addr;
  *__prior_counter_ptr = 0;
  mmap_addr += sizeof(long long);
  __current_counter_ptr = (volatile long long *) mmap_addr;
  *__current_counter_ptr = 0;
  mmap_addr += sizeof(long long);           
#endif
}
#endif  
