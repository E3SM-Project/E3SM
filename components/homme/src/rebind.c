#ifdef _AIX


#include <sys/processor.h>
#include <sys/thread.h>
#include <stdio.h>

static cpu_t my_cpu;
static tid_t mytid;

int unbind()
{
  int rc;

  my_cpu = mycpu();
  mytid = thread_self ();

  rc = bindprocessor(BINDTHREAD, mytid, PROCESSOR_CLASS_ANY);
  if(rc != 0) 
    fprintf(stderr,"Failed to unbind thread %d from cpu %d\n",mytid,mycpu);

  
}

int rebind()
{
  int rc;
  rc = bindprocessor(BINDTHREAD, mytid, my_cpu);
  if(rc != 0) 
    fprintf(stderr,"Failed to bind thread %d to cpu %d\n",mytid,mycpu);
}


#endif
