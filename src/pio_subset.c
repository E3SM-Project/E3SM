#include <pio.h>
#include <pio_internal.h>

typedef struct mapsort
{
  int rfrom;
  int soffset;
  PIO_Offset iomap;
} mapsort;



int compare_offsets(const void *a,const void *b) 
{
  mapsort *x = (mapsort *) a;
  mapsort *y = (mapsort *) b;
  return (int) (x->iomap - y->iomap);
}    



int subset_rearrange_create(iosystem_desc_t *ios,const int maplen, const PIO_Offset compmap[], 
			    const int gsize[], const int ndims, const int nioproc, io_desc_t *iodesc)
{
  int *dest_ioproc;
  PIO_Offset *dest_ioindex;
  int tsize;
  MPI_Datatype dtype;
  int taskratio;
  int nprocs = ios->num_comptasks;
  int nioprocs = ios->num_iotasks;
  int *sndlths; 
  int *sdispls;
  int *recvlths;
  int *rdispls;
  MPI_Datatype *dtypes;
  int i, j, jlast;
  bool hs=false;
  bool isend=false;
  PIO_Offset *iomap=NULL;
  int ierr = PIO_NOERR;
  int iomaplen=0;
  mapsort *map;
  PIO_Offset *destoffset=NULL;


  printf("maplen %ld nioproc %d\n",maplen, nioproc); 


  dest_ioproc = (int *) malloc(maplen*sizeof(int));
  dest_ioindex = (PIO_Offset *) malloc(maplen*sizeof(PIO_Offset));

  iodesc->ndof = maplen;

  PIO_Offset_size(&dtype, &tsize);
  
  taskratio = nprocs/nioprocs;

  // Each compute task sends to only one IO task
  for(i=0;i<maplen;i++){
    dest_ioproc[i] = ios->comp_rank/taskratio;
  }

  // Pass the maplen from each compute task to its associated IO task
  
  sndlths = (int *) malloc(nprocs*sizeof(int)); 
  sdispls= (int *) malloc(nprocs*sizeof(int));
  recvlths= (int *) malloc(nprocs*sizeof(int));
  rdispls= (int *) malloc(nprocs*sizeof(int));
  dtypes= (MPI_Datatype *) malloc(nprocs*sizeof(MPI_Datatype));

  for(i=0;i<nprocs;i++){
    sndlths[i] = 0;
    sdispls[i] = 0;
    recvlths[i] = 0;
    rdispls[i] = 0;
    dtypes[i] = MPI_INT;
  }
  sndlths[ ios->ioranks[ dest_ioproc[0] ] ] = 1;
  if(ios->ioproc){
    for(i=0;i<nprocs;i++)
      if(ios->io_rank == i/taskratio){
	recvlths[i]=1;
	rdispls[i] = sizeof(int) * (i % taskratio);
      }
  }
  int recvlen[taskratio];
  int jstart[taskratio];

  pio_swapm(&maplen, sndlths, sdispls, dtypes, 
	    recvlen, recvlths, rdispls, dtypes, 
	    ios->union_comm, hs, isend, MAX_GATHER_BLOCK_SIZE);


  // Now pass the map
  if(ios->ioproc){
    for(i=0;i<taskratio;i++){
      iomaplen+=recvlen[i];
      printf("%d recvlen=%d\n",i,recvlen[i]);
    }
    map = (mapsort *) malloc(iomaplen * sizeof(mapsort));    
    iomap = (PIO_Offset *) calloc(iomaplen,sizeof(PIO_Offset));
    destoffset = (PIO_Offset *) calloc(iomaplen,sizeof(PIO_Offset));
  }
  sndlths[ ios->ioranks[ dest_ioproc[0] ] ] = maplen;
  if(ios->ioproc){
    jlast =0;
    for(i=0;i<nprocs;i++)
      if(ios->io_rank == i/taskratio){
	recvlths[i]= recvlen[ i%taskratio];
	jstart[ i%taskratio ] = jlast;
	printf("here %d %d jstart %d\n",i,i%taskratio, jlast);
	for(j=jlast;j<jlast+recvlths[i];j++){
	  (map+j)->rfrom = i;
	  (map+j)->soffset = j-jlast;
	}
	
	jlast = j;
      }
    
  }
  for(i=0;i<nprocs;i++)
    dtypes[i] = dtype;

  pio_swapm(compmap, sndlths, sdispls, dtypes, 
	    iomap, recvlths, rdispls, dtypes, 
	    ios->union_comm, hs, isend, MAX_GATHER_BLOCK_SIZE);

  if(ios->ioproc){
    for(i=0;i<iomaplen;i++){
      (map+i)->iomap = iomap[i];
    }
    for(i=0;i<iomaplen;i++)
      printf("before from %d offset %d iomap %ld\n",(map+i)->rfrom,(map+i)->soffset,(map+i)->iomap);

    // sort the mapping, this will transpose the data into IO order        
    qsort(map, iomaplen, sizeof(mapsort), compare_offsets); 

    // Now we want to send the iotask location back to the compute task
    int destloc;
    for(i=0;i<iomaplen;i++){
      iomap[i]=(map+i)->iomap;  
      destloc = jstart[(map+i)->rfrom % taskratio] + ((map+i)->soffset)  ;
      if(destloc<0 || destloc>iomaplen)
	fprintf(stderr,"Error computing destloc %d \n",destloc);

      destoffset[ destloc ] = i;
      printf("%d %d %ld %d %d\n",i,destloc, destoffset[destloc], (map+i)->rfrom, (map+i)->soffset);
      //      printf("after from %d offset %d iomap %ld \n",(map+i)->rfrom,(map+i)->soffset,(map+i)->iomap);
    }      



  }
  pio_swapm(destoffset, recvlths, rdispls, dtypes, 
	    dest_ioindex, sndlths, sdispls, dtypes, 
	    ios->union_comm, hs, isend, MAX_GATHER_BLOCK_SIZE);
  for(i=0;i<maplen;i++)
    printf("i=%d dest_ioindex = %ld\n",i,dest_ioindex[i]);






  return ierr;


}
  
    
  
