#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#if defined(BGL) || defined(BGP)

#include "mpi.h"
#include <bglpersonality.h>
#include "rts.h"
#include "bginfo.h"

int rank;
int np;
int my_name_len;
char my_name[255];

void identity(MPI_Fint *comm, int *iotask){

   
   MPI_Comm comm2;
   comm2 = MPI_Comm_f2c(*comm);
   MPI_Comm_rank(comm2,&rank);
   MPI_Comm_size(comm2,&np);
   MPI_Get_processor_name(my_name, &my_name_len);

   /*  Get the personality  */
   BGLPersonality pers;
   char message[100];

   rts_get_personality(&pers, sizeof(pers));
   BGLPersonality_getLocationString(&pers, message);
    
   int numIONodes,numPsets,numNodesInPset,rankInPset;
   numIONodes = BGLPersonality_numIONodes(&pers);
   numPsets = BGLPersonality_numPsets(&pers);
   numNodesInPset = BGLPersonality_numNodesInPset(&pers);
   rankInPset = BGLPersonality_rankInPset(&pers);
    
/*
   if(rank == 0) { printf("number of IO nodes in block: %i \n",numIONodes);}
   if(rank == 0) { printf("number of Psets in block : %i \n",numPsets);}
   if(rank == 0) { printf("number of compute nodes in Pset: %i \n",numNodesInPset);}
   printf("MPI task %i is rank %i in Pset \n",rank, rankInPset);
*/

   int psetNum;
   psetNum = BGLPersonality_psetNum(&pers);
   if((*iotask)>0) {
      printf( "%04i (%-50s %s) %i yes\n", rank, my_name, message, psetNum );
   } else {
      printf( "%04i (%-50s %s) %i --\n", rank, my_name, message, psetNum);
   }

}

void determineiotasks(MPI_Fint *comm, int *numiotasks,int *base, int *stride, int *rearr, int *iamIOtask){

/*  Some concepts:

     processor set: 	A group of processors on the Blue Gene system which have 
			one or more IO processor (Pset)

     IO-node:		A special Blue Gene node dedicated to performing IO.  There 
			are one or more per processor set

     IO-client:		This is software concept.  This refers to the MPI task 
			which performs IO for the PIO library 
*/
   int psetNum;                                 
   int coreId;
   int iam;
   MPI_Comm comm2;
   comm2 = MPI_Comm_f2c(*comm);
   MPI_Comm_rank(comm2,&rank);
   MPI_Comm_size(comm2,&np);
   MPI_Get_processor_name(my_name, &my_name_len);

   /*  Get the personality  */
   BGLPersonality pers;
   char message[100];

if((*rearr) > 0) {
   rts_get_personality(&pers, sizeof(pers));

   int numIONodes,numPsets,numNodesInPset,rankInPset;
   int numiotasks_per_node,remainder,numIONodes_per_pset;
   int lstride;

   /* total number of IO-nodes */
   numIONodes = BGLPersonality_numIONodes(&pers);

   /* Number of computational nodes in processor set */
   numNodesInPset = BGLPersonality_numNodesInPset(&pers);

   
   if((*numiotasks) < 0 ) { 
       /* negative numiotasks value indicates that this is the number per IO-node */
       (*numiotasks) = - (*numiotasks);
       if((*numiotasks) > numNodesInPset) {
          numiotasks_per_node = numNodesInPset;
       } else  {
          numiotasks_per_node = (*numiotasks);
       }
       remainder = 0;
   } else if ((*numiotasks) > 0 ) {
       /* balance the number of iotasks to number of IO nodes */
       numiotasks_per_node = floor((float)(*numiotasks)/ (float) numIONodes);
       remainder = (*numiotasks) - numiotasks_per_node * numIONodes;
   } else if ((*numiotasks) == 0 ) {
       if((*stride) > 0) {
          numiotasks_per_node = numNodesInPset/(*stride);
       } else {
          numiotasks_per_node = 6;  /* default number of IO-client per IO-node is not otherwise specificied */
       }
       remainder = 0;
   } 

   /* number of IO nodes with a larger number of io-client per io-node */
   if(remainder > 0) {
       if(rank ==0) {printf("Unbalanced IO-configuration: %i IO-nodes have %i IO-clients : %i IO-nodes have %i IO-clients \n",
		remainder, numiotasks_per_node+1, numIONodes-remainder,numiotasks_per_node);}
       lstride = floor((float)numNodesInPset/(float)(numiotasks_per_node+1));
   } else {
       if(rank == 0) {printf("Balanced IO-configuration: %i IO-nodes have %i IO-clients\n",numIONodes-remainder, numiotasks_per_node);}
       lstride = floor((float)numNodesInPset/(float)numiotasks_per_node);
   }
  
   /* Number of processor sets */
   numPsets = BGLPersonality_numPsets(&pers);
 
   /* number of IO nodes in processor set (I need to add
      code to deal with the case where numIONodes_per_pset != 1 works 
      correctly) */
   numIONodes_per_pset = numIONodes/numPsets;

  /* Determine which core on node....  I don't want to put more than one io-task per node */
   coreId = rts_get_processor_id();

   /* What is the rank of this node in the processor set */
   rankInPset = BGLPersonality_rankInPset(&pers);

   /* determine the processor set that this node belongs to */
   psetNum = BGLPersonality_psetNum(&pers);

/* printf("Pset #: %i has %i nodes in Pset\n",psetNum,numNodesInPset); */

   (*iamIOtask) = 0;   /* initialize to zero */

   if((*stride) == 1) (*base) = 0;  /* Reset the base to 0 if we are using all tasks */
   /* start stridding MPI tasks from base task */ 
   iam = rankInPset-(*base);
   if (iam >= 0)  {
       /* mark tasks that will be IO-tasks  or IO-clients */
       if((iam % lstride == 0) && (coreId == 0) ) {  /* only io tasks indicated by stride and coreId = 0 */
           if((iam/lstride) < numiotasks_per_node) { 
	      /* only set the first (numiotasks_per_node - 1) tasks */
              (*iamIOtask) = 1;
           } else if ((iam/lstride) == numiotasks_per_node) {
	      /*  If there is an uneven number of io-clients to io-nodes 
	          allocate the first remainder - 1 processor sets to 
	          have a total of numiotasks_per_node */
	      if(psetNum < remainder) {(*iamIOtask) = 1;};   
           }
       }
   }
   (*numiotasks) = numiotasks_per_node * numIONodes;
}  else {
   /* We are not doing rearrangement.... so all tasks are io-tasks */
   (*iamIOtask) = 1;
}

}

#endif
