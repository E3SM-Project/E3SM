#if defined (IRIX64) && defined(PIN_CPUS)
#define SN0 1
#define SN0XXL 1
#define _KMEMUSER 1

#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/syssgi.h>
#include <sys/pmo.h>
#include <sys/sysmp.h>
#include <sys/nodemask.h>
#include <sys/SN/SN0/arch.h>
#include <alloca.h>
#include <strings.h>
#include <errno.h>

/* #include <sys/miser_public.h> */
#include <cpuset.h>


/**********************************************************************/
/* cpuset info */

id_type_t
__get_current_cpuset_name(void)
{
  cpuset_request_t req;
  /* cpuset_qname_t *names = (cpuset_qname_t*) &(req.csr_u.cs_qname); */
 
  req.request = CPUSET_QUERY_CURRENT;
  if (sysmp(MP_CPUSET, &req) == -1) {
    /* Not in a cpuset */
    return -1;
  }
  /* return  (id_type_t) (names->qname); */
  return  (id_type_t) (req.csr_u.cs_qname.qname);
}



int
__get_cpus_in_cpuset(id_type_t cpuset_id, cpuid_t *cpuset)
{
  cpuset_request_t req;
  cpuset_queue_t* cs = (cpuset_queue_t*) &(req.csr_u.cs_queue);

  cs->queue = cpuset_id;
  cs->cpuid = (uint64_t) cpuset;

  req.request = CPUSET_QUERY_CPUS;
  if (sysmp(MP_CPUSET, &req)) {
    fprintf(stderr, "Could not get status for cpuset 0x%llx\n",
                 (unsigned long long int) cpuset_id);
    perror("sysmp");
    return -1;
  }

  return cs->count;
}



/* cpuset should already be zeroed on input
** (since we don't know how long it is, we don't do it)
*/
int
__get_cpus_in_my_cpuset(cpuid_t *cpuset)
{
  id_type_t cpuset_name;
  int count = 0;

  cpuset_name = __get_current_cpuset_name();
  if (cpuset_name != (id_type_t) -1) {
    count = __get_cpus_in_cpuset(cpuset_name, cpuset);
    if (count < 0) count = 0;
  }
  return count;
}

  


/*************************************************************************/
/* nodemask info */

int
__get_nodes_in_my_nodemask(char *nodemask)
{
  cnodemask_t sys_nodemask;
  int i, rtn_value = 0;

  bzero(&sys_nodemask, sizeof(sys_nodemask));

  if(pmoctl(PMO_GETNODEMASK_UINT64,&sys_nodemask,sizeof(sys_nodemask)) < 0) {
    perror("pmoctl(PMO_GETNODEMASK_UINT64)");
    exit(1);
  }

  for (i = 0; i < MAX_COMPACT_NODES; i++) {
    if (CNODEMASK_TSTB(sys_nodemask, i) != 0) {
      nodemask[i] = 1;
      rtn_value++;
    }
  }

  return rtn_value;
}


/**************************************************************************/





void
__compute_placement(int my_rank, int *cpu, int *memory)
{
  long ncpus = sysmp(MP_NPROCS);
  char *nodemask = alloca((unsigned int)ncpus);
  cnodeid_t *cpu_node_mapping = alloca((unsigned int)(ncpus*sizeof(cnodeid_t)));
  cpuid_t  *cpuset = alloca((unsigned int)(ncpus*sizeof(cpuid_t)));
  int i, cpuset_size;

  *cpu = -1;
  *memory = -1;
  bzero(nodemask, ncpus);
  bzero(cpu_node_mapping, ncpus*sizeof(cnodeid_t));


  /* Construct a list of the cpus assigned to this job.  If we
  ** are in a miser_cpuset, use that.  If not, use the nodemask
  ** to construct the list (all jobs always have a nodemask).
  **
  **   Once constructed, simply index the list with my_rank.
  ** We index from the top down, so that if the cpu list is larger
  ** than the number of cpus we are actually using, the unused cpus
  ** will be the smaller numbered ones.  This helps the default case
  ** by saving cpu-0 for last (and hopefully, not using it at all).
  */

  /* Ask the kernel for the cpu -> node mapping */
  if (sysmp(MP_NUMA_GETCPUNODEMAP,
              (void *)cpu_node_mapping, sizeof(cnodeid_t) * ncpus) != 0)
  {
    perror("Could not get cpu->node mapping sysmp(MP_NUMA_GETCPUNODEMAP)");
    exit(1);
  }


  if ((cpuset_size = __get_cpus_in_my_cpuset(cpuset)) <= 0) {
    /* not in a cpuset */
    if (__get_nodes_in_my_nodemask(nodemask) <= 0) {
      fprintf(stderr, "Warning: invalid nodemask in __compute_placement\n");
      return;
    }
    else {
      /* Make a compact list of cpus under the nodemask */
      for (i = 0, cpuset_size = 0; i < ncpus; i++) {
        if (nodemask[cpu_node_mapping[i]]) {
          cpuset[cpuset_size] = i;
          cpuset_size++;
        }
      }
      if (cpuset_size <= 0) {
        fprintf(stderr, "Warning: No cpus found under nodemask ???\n");
        return;
      }
    }
  }


  if (my_rank >= cpuset_size) {  
    fprintf(stderr, "Warning: more mp processes than available cpus (%d)\n", 
                cpuset_size);
    my_rank = my_rank % cpuset_size;
  }  

  *cpu = cpuset[my_rank];
  *memory = cpu_node_mapping[*cpu];

}




/* Note that this is a general routine; the cpu and the memory  */
/* do NOT need to be on the same node (although that is faster) */
void
__place_process(int cpu_to_use, int memory_to_use)
{
  pmo_handle_t mld, mldset;
  raff_info_t affinity_info;
  char name[PATH_MAX+1];
  pid_t my_pid = getpid();
  int status;

  /* First, create an mld (a singleton mldset) ... */
  if((mld = mld_create(0, 16*1024)) < 0) {
    perror("mld_create(0, 16*1024)");
    exit(1);
  }
  if ((mldset = mldset_create(&mld, 1)) < 0) {
    perror("mldset_create(&mld, 1)");
    exit(1);
  }

  /* ... now place the mld onto the requested node ... */
  sprintf(name, "/hw/nodenum/%d", memory_to_use);
  affinity_info.resource = name;
  affinity_info.reslen = (unsigned short) strlen(name);
  affinity_info.restype = RAFFIDT_NAME;
  affinity_info.radius = 0;
  affinity_info.attr = RAFFATTR_ATTRACTION;

  status = mldset_place(mldset, TOPOLOGY_PHYSNODES,
			&affinity_info, 1, RQMODE_MANDATORY);
  if (status != 0) {
    /* If MANDATORY fails, try again with ADVISORY */
    status = mldset_place(mldset, TOPOLOGY_PHYSNODES,
			&affinity_info, 1, RQMODE_ADVISORY);
    if (status != 0) {
      perror("mldset_place(RQMODE_ADVISORY)");
      exit(1);
    }
  }

  /* ... and link the process to the mld. */
  if (process_mldlink(my_pid, mld, RQMODE_MANDATORY) < 0) {
    if (process_mldlink(my_pid, mld, RQMODE_ADVISORY) < 0) {
      perror("process_mldlink(RQMODE_ADVISORY)");
      exit(1);
    }
  }

  /* Now force the process to run on the requested cpu */
  sysmp(MP_RUNANYWHERE);  /* Break any previous affinity */
  if (sysmp(MP_MUSTRUN, cpu_to_use) < 0) {
    perror("sysmp(MP_MUSTRUN)");
    exit(1);
  }

}




void
__mp_assign_to_cpu(uint32_t relative_cpu)
{
  int cpu_to_use, memory_to_use;

  __compute_placement(relative_cpu, &cpu_to_use, &memory_to_use);

  if ((cpu_to_use == -1) || (memory_to_use == -1)) return;

  __place_process(cpu_to_use, memory_to_use);
}


/* The FORTRAN entry point */
void
mp_assign_to_cpu_(uint32_t *relative_cpu)
{ __mp_assign_to_cpu(*relative_cpu); }



void
oinker_(void)
{
  sysmp(MP_RUNANYWHERE);  /* Break any previous affinity */
}


void
unoinker_(void)
{
  int cpu_to_use, memory_to_use;
  __compute_placement(0, &cpu_to_use, &memory_to_use);

  /* Now force the process to run on the requested cpu */
  sysmp(MP_RUNANYWHERE);  /* Break any previous affinity */
  if (sysmp(MP_MUSTRUN, cpu_to_use) < 0) {
    perror("sysmp(MP_MUSTRUN)");
    exit(1);
  }
}
#endif
