#include <pio.h>
#include <pio_internal.h>

PIO_Offset PIO_BUFFER_SIZE_LIMIT= 100000000; // 100MB default limit

 // Changes to PIO_BUFFER_SIZE_LIMIT only apply to files opened after the change
 PIO_Offset PIOc_set_buffer_size_limit(const PIO_Offset limit)
 {
   PIO_Offset oldsize; 
   oldsize = PIO_BUFFER_SIZE_LIMIT;
   if(limit>0)
     PIO_BUFFER_SIZE_LIMIT=limit;
   return(oldsize);
 }

 int pio_write_darray_nc(file_desc_t *file, io_desc_t *iodesc, const int vid, void *IOBUF, void *fillvalue)
 {
   iosystem_desc_t *ios;
   var_desc_t *vdesc;
   int ndims;
   int ierr;
   int i;
   int msg;
   int mpierr;
   int dsize;
   MPI_Status status;
   PIO_Offset usage;
   int request;
   int fndims;
   PIO_Offset tdsize;

   tdsize=0;
   ierr = PIO_NOERR;

   ios = file->iosystem;
   if(ios == NULL){
     fprintf(stderr,"Failed to find iosystem handle \n");
     return PIO_EBADID;
   }
   vdesc = (file->varlist)+vid;

   if(vdesc == NULL){
     fprintf(stderr,"Failed to find variable handle %d\n",vid);
     return PIO_EBADID;
   }
   ndims = iodesc->ndims;
   msg = 0;

   if(ios->async_interface && ! ios->ioproc){
     if(ios->comp_rank==0) 
       mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
     mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
   }

   ierr = PIOc_inq_varndims(file->fh, vid, &fndims);

   if(ios->ioproc){
     io_region *region;
     int ncid = file->fh;
     int regioncnt;
     int rrcnt;
     void *bufptr;
     void *tmp_buf=NULL;
     int tsize;
     size_t start[fndims];
     size_t count[fndims];
     int buflen, j;

     PIO_Offset *startlist[iodesc->maxregions];
     PIO_Offset *countlist[iodesc->maxregions];


#ifdef _MPISERIAL
     tsize = iodesc->basetype;
#else
     MPI_Type_size(iodesc->basetype, &tsize);
#endif
     region = iodesc->firstregion;

     if(vdesc->record >= 0 && ndims<fndims)
       ndims++;
#ifdef _PNETCDF
     if(file->iotype == PIO_IOTYPE_PNETCDF){
       // make sure we have room in the buffer ;
	 ierr = ncmpi_inq_buffer_usage(ncid, &usage);
	 usage += tsize*(iodesc->maxiobuflen);
	 MPI_Allreduce(MPI_IN_PLACE, &usage, 1,  MPI_LONG_LONG,  MPI_MAX, ios->io_comm);
	 //         if(ios->io_rank==0) printf("%s %d %d %d\n",__FILE__,__LINE__,iodesc->maxregions,usage);
	 if(usage >= PIO_BUFFER_SIZE_LIMIT){
	   flush_output_buffer(file);
	 }

     }
#endif

     rrcnt=0;
     for(regioncnt=0;regioncnt<iodesc->maxregions;regioncnt++){
       for(i=0;i<ndims;i++){
	 start[i] = 0;
	 count[i] = 0;
       }
       if(region != NULL){
	 bufptr = (void *)((char *) IOBUF+tsize*region->loffset);
	 // this is a record based multidimensional array
	 if(vdesc->record >= 0){
	   start[0] = vdesc->record;
	   for(i=1;i<ndims;i++){
	     start[i] = region->start[i-1];
	     count[i] = region->count[i-1];
	   }
	  if(count[1]>0)
	    count[0] = 1;
	 // Non-time dependent array
	 }else{
	   for( i=0;i<ndims;i++){
	     start[i] = region->start[i];
	     count[i] = region->count[i];
	   }
	 }
       }

       switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
       case PIO_IOTYPE_NETCDF4P:
	 ierr = nc_var_par_access(ncid, vid, NC_COLLECTIVE);
	 switch(iodesc->basetype){
	 case MPI_DOUBLE:
	 case MPI_REAL8:
	   ierr = nc_put_vara_double (ncid, vid,(size_t *) start,(size_t *) count, (const double *) bufptr); 
	   break;
	 case MPI_INTEGER:
	   ierr = nc_put_vara_int (ncid, vid, (size_t *) start, (size_t *) count, (const int *) bufptr); 
	   break;
	 case MPI_FLOAT:
	 case MPI_REAL4:
	   ierr = nc_put_vara_float (ncid, vid, (size_t *) start, (size_t *) count, (const float *) bufptr); 
	   break;
	 default:
	   fprintf(stderr,"Type not recognized %d in pioc_write_darray\n",(int) iodesc->basetype);
	 }
	 break;
       case PIO_IOTYPE_NETCDF4C:
#endif
       case PIO_IOTYPE_NETCDF:
	 {
#ifdef _MPISERIAL
	   dsize = iodesc->basetype;
#else
	   mpierr = MPI_Type_size(iodesc->basetype, &dsize);
 #endif
	   size_t tstart[ndims], tcount[ndims];
	   if(ios->io_rank==0){
         // FIX(SPM, 100714)  don't use i, use something like myrank and iam
	     for(i=0;i<iodesc->num_aiotasks;i++){
	       if(i==0){	    
		 buflen=1;
		 for(j=0;j<ndims;j++){
		   tstart[j] =  start[j];
		   tcount[j] =  count[j];
		   buflen *= tcount[j];
		   tmp_buf = bufptr;
		 }
	       }else{
		 mpierr = MPI_Send( &ierr, 1, MPI_INT, i, 0, ios->io_comm);  // handshake - tell the sending task I'm ready
		 mpierr = MPI_Recv( &buflen, 1, MPI_INT, i, 1, ios->io_comm, &status);
		 if(buflen>0){
		   mpierr = MPI_Recv( tstart, ndims, MPI_OFFSET, i, ios->num_iotasks+i, ios->io_comm, &status);
		   mpierr = MPI_Recv( tcount, ndims, MPI_OFFSET, i,2*ios->num_iotasks+i, ios->io_comm, &status);
		   tmp_buf = malloc(buflen * dsize);	
		   mpierr = MPI_Recv( tmp_buf, buflen, iodesc->basetype, i, i, ios->io_comm, &status);
		 }
	       }

	       if(buflen>0){
		 if(iodesc->basetype == MPI_INTEGER){
		   ierr = nc_put_vara_int (ncid, vid, tstart, tcount, (const int *) tmp_buf); 
		 }else if(iodesc->basetype == MPI_DOUBLE || iodesc->basetype == MPI_REAL8){
		   ierr = nc_put_vara_double (ncid, vid, tstart, tcount, (const double *) tmp_buf); 
		 }else if(iodesc->basetype == MPI_FLOAT || iodesc->basetype == MPI_REAL4){
		   ierr = nc_put_vara_float (ncid,vid, tstart, tcount, (const float *) tmp_buf); 
		 }else{
		   fprintf(stderr,"Type not recognized %d in pioc_write_darray\n",(int) iodesc->basetype);
		 }
		 if(ierr == PIO_EEDGE){
		   for(i=0;i<ndims;i++)
		     fprintf(stderr,"dim %d start %ld count %ld\n",i,tstart[i],tcount[i]);
		 }
		 if(tmp_buf != bufptr)
		   free(tmp_buf);
	       }
	     }
	   }else if(ios->io_rank < iodesc->num_aiotasks ){
	     buflen=1;
	     for(i=0;i<ndims;i++){
	       tstart[i] = (size_t) start[i];
	       tcount[i] = (size_t) count[i];
	       buflen*=tcount[i];
	       //               printf("%s %d %d %d %d\n",__FILE__,__LINE__,i,tstart[i],tcount[i]);
	     }
	     //	     printf("%s %d %d %d %d %d %d %d %d %d\n",__FILE__,__LINE__,ios->io_rank,tstart[0],tstart[1],tcount[0],tcount[1],buflen,ndims,fndims);
	     mpierr = MPI_Recv( &ierr, 1, MPI_INT, 0, 0, ios->io_comm, &status);  // task0 is ready to recieve
	     mpierr = MPI_Rsend( &buflen, 1, MPI_INT, 0, 1, ios->io_comm);
	     if(buflen>0) {
	       mpierr = MPI_Rsend( tstart, ndims, MPI_OFFSET, 0, ios->num_iotasks+ios->io_rank, ios->io_comm);
	       mpierr = MPI_Rsend( tcount, ndims, MPI_OFFSET, 0,2*ios->num_iotasks+ios->io_rank, ios->io_comm);
	       mpierr = MPI_Rsend( bufptr, buflen, iodesc->basetype, 0, ios->io_rank, ios->io_comm);
	     }
	   }
	   break;
	 }
	 break;
 #endif
 #ifdef _PNETCDF
       case PIO_IOTYPE_PNETCDF:
	 for( i=0,dsize=1;i<ndims;i++){
	   dsize*=count[i];
	 }
	 tdsize += dsize;
	 //	 if(dsize==1 && ndims==2)
	 //	 printf("%s %d %d\n",__FILE__,__LINE__,iodesc->basetype);

	 /*	 if(regioncnt==0){
	   for(i=0;i<iodesc->maxregions;i++){
	     startlist[i] = (PIO_Offset *) calloc(fndims, sizeof(PIO_Offset));
	     countlist[i] = (PIO_Offset *) calloc(fndims, sizeof(PIO_Offset));
	   }
	 }
	 */
	 if(dsize>0){
	   //	   printf("%s %d %d %d\n",__FILE__,__LINE__,ios->io_rank,dsize);
	   startlist[rrcnt] = (PIO_Offset *) calloc(fndims, sizeof(PIO_Offset));
	   countlist[rrcnt] = (PIO_Offset *) calloc(fndims, sizeof(PIO_Offset));
	   for( i=0; i<fndims;i++){
	     startlist[rrcnt][i]=start[i];
	     countlist[rrcnt][i]=count[i];
	   }
	   rrcnt++;
	 }	 
	 if(regioncnt==iodesc->maxregions-1){
	   //printf("%s %d %d %ld %ld\n",__FILE__,__LINE__,ios->io_rank,iodesc->llen, tdsize);
	   //	   ierr = ncmpi_put_varn_all(ncid, vid, iodesc->maxregions, startlist, countlist, 
	   //			     IOBUF, iodesc->llen, iodesc->basetype);
	   
	   //printf("%s %d %ld \n",__FILE__,__LINE__,IOBUF);
	   ierr = ncmpi_bput_varn(ncid, vid, rrcnt, startlist, countlist, 
				  IOBUF, iodesc->llen, iodesc->basetype, &request);
	   pio_push_request(file,request);
	   for(i=0;i<rrcnt;i++){
	     free(startlist[i]);
	     free(countlist[i]);
	   }
	 }
	 break;
#endif
       default:
	 ierr = iotype_error(file->iotype,__FILE__,__LINE__);
       }
       if(region != NULL)
	 region = region->next;
     } //    for(regioncnt=0;regioncnt<iodesc->maxregions;regioncnt++){
   } // if(ios->ioproc)

   ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

   return ierr;
 }

int pio_write_darray_multi_nc(file_desc_t *file, io_desc_t *iodesc,const int nvars, const int vid[], 
			      void *IOBUF, const int frame[], void *fillvalue[])
 {
   iosystem_desc_t *ios;
   var_desc_t *vdesc;
   int ndims;
   int ierr;
   int i;
   int msg;
   int mpierr;
   int dsize;
   MPI_Status status;
   PIO_Offset usage;
   int request;
   int fndims;
   PIO_Offset tdsize;
   int tsize;
   int ncid;
   tdsize=0;
   ierr = PIO_NOERR;

   ios = file->iosystem;
   if(ios == NULL){
     fprintf(stderr,"Failed to find iosystem handle \n");
     return PIO_EBADID;
   }
   vdesc = (file->varlist)+vid[0];
   ncid = file->fh;

   if(vdesc == NULL){
     fprintf(stderr,"Failed to find variable handle %d\n",vid[0]);
     return PIO_EBADID;
   }
   ndims = iodesc->ndims;
   msg = 0;

   if(ios->async_interface && ! ios->ioproc){
     if(ios->comp_rank==0) 
       mpierr = MPI_Send(&msg, 1,MPI_INT, ios->ioroot, 1, ios->union_comm);
     mpierr = MPI_Bcast(&(file->fh),1, MPI_INT, ios->compmaster, ios->intercomm);
   }

   ierr = PIOc_inq_varndims(file->fh, vid[0], &fndims);
#ifdef _MPISERIAL
   tsize = iodesc->basetype;
#else
   MPI_Type_size(iodesc->basetype, &tsize);
#endif



   if(iodesc->needsfill){
     PIO_Offset start[fndims];
     PIO_Offset count[fndims];
     void *bufptr=NULL;
     size_t totalsize=0;

     printf("%s %d \n",__FILE__,__LINE__);

     if(ios->io_rank==0){
       totalsize=1;
       for(int id=(fndims-ndims);id<fndims;id++){
	 start[id]=0;
	 count[id]=iodesc->gsize[id-(fndims-ndims)];
	 totalsize*=count[id];
       }
       bufptr = malloc(tsize*totalsize);
     }else{
       for(int id=0;id<fndims;id++){
	 start[id]=0;
	 count[id]=0;
       }
     }

     for(int nv=0;nv<nvars;nv++){
       printf("%s %d %d\n",__FILE__,__LINE__,nv);
       vdesc = (file->varlist)+vid[nv];
       if(fndims>ndims && ios->io_rank==0){
	 start[0]=vdesc->record;
	 count[0]=1;
       }
       switch(iodesc->basetype){
       case MPI_DOUBLE:
       case MPI_REAL8:
	 for(i=0;i<totalsize;i++){
	   ((double *)bufptr)[i]  = ((double *) fillvalue)[nv];
	 }
	 PIOc_put_vara_double(ncid, vid[nv], start, count, (const double *) bufptr);
	 break;
       case MPI_INTEGER:
         printf("%s %d %d %d %d\n",__FILE__,__LINE__,vid[nv],((int *) fillvalue)[nv],totalsize);
	 for(i=0;i<totalsize;i++){
	   ((int *)bufptr)[i]  = ((int *) fillvalue)[nv];
	 }
	 PIOc_put_vara_int(ncid, vid[nv], start, count, (const int *) bufptr);
	 break;
       case MPI_FLOAT:
       case MPI_REAL4:
	 for(i=0;i<totalsize;i++){
	   ((float *)bufptr)[i]  = ((float *) fillvalue)[nv];
	 }
	 PIOc_put_vara_float(ncid, vid[nv], start, count, (const float *) bufptr);
	 break;
       default:
	 fprintf(stderr,"Type not recognized %d in pioc_write_darray\n",(int) iodesc->basetype);
       }	 
     }
#ifdef _PNETCDF
     // Need to flush the buffer of fill values before the next write
     printf("%s %d \n",__FILE__,__LINE__);
     if(ios->ioproc && file->iotype == PIO_IOTYPE_PNETCDF){
       flush_output_buffer(file);
     }
#endif
     if(totalsize>0){
       free(bufptr);
     }
   }



   if(ios->ioproc){
     io_region *region;
     int regioncnt;
     int rrcnt;
     void *bufptr;
      int buflen, j;
     size_t start[fndims];
     size_t count[fndims];


     PIO_Offset *startlist[iodesc->maxregions];
     PIO_Offset *countlist[iodesc->maxregions];

     ncid = file->fh;
     region = iodesc->firstregion;

     if(vdesc->record >= 0 && ndims<fndims)
       ndims++;
#ifdef _PNETCDF
     if(file->iotype == PIO_IOTYPE_PNETCDF){
       // make sure we have room in the buffer ;
	 ierr = ncmpi_inq_buffer_usage(ncid, &usage);
	 usage += nvars*tsize*(iodesc->maxiobuflen);
	 printf("%s %d %ld %ld\n",__FILE__,__LINE__,usage,PIO_BUFFER_SIZE_LIMIT);
	 MPI_Allreduce(MPI_IN_PLACE, &usage, 1,  MPI_LONG_LONG,  MPI_MAX, ios->io_comm);
	 //	 printf("%s %d %ld %ld %ld\n",__FILE__,__LINE__,usage,PIO_BUFFER_SIZE_LIMIT,nvars*tsize*(iodesc->maxiobuflen));
	 if(usage >= PIO_BUFFER_SIZE_LIMIT){	   
	   flush_output_buffer(file);
	 }

     }
#endif

     rrcnt=0;
     for(regioncnt=0;regioncnt<iodesc->maxregions;regioncnt++){
       for(i=0;i<ndims;i++){
	 start[i] = 0;
	 count[i] = 0;
       }
       if(region != NULL){
	 // this is a record based multidimensional array
	 if(vdesc->record >= 0){
	   start[0] = frame[0];
	   for(i=1;i<ndims;i++){
	     start[i] = region->start[i-1];
	     count[i] = region->count[i-1];
	   }
	  if(count[1]>0)
	    count[0] = 1;
	 // Non-time dependent array
	 }else{
	   for( i=0;i<ndims;i++){
	     start[i] = region->start[i];
	     count[i] = region->count[i];
	   }
	 }
       }

       switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
       case PIO_IOTYPE_NETCDF4P:
	 for(int nv=0; nv<nvars; nv++){
	   if(vdesc->record >= 0){
	     start[0] = frame[nv];
	   }
	   if(region != NULL){
	     bufptr = (void *)((char *) IOBUF + tsize*(nv*iodesc->llen + region->loffset));
	   }
	   ierr = nc_var_par_access(ncid, vid[nv], NC_COLLECTIVE);

	   switch(iodesc->basetype){
	   case MPI_DOUBLE:
	   case MPI_REAL8:
	     ierr = nc_put_vara_double (ncid, vid[nv],(size_t *) start,(size_t *) count, (const double *) bufptr); 
	     break;
	   case MPI_INTEGER:
	     ierr = nc_put_vara_int (ncid, vid[nv], (size_t *) start, (size_t *) count, (const int *) bufptr); 
	     break;
	   case MPI_FLOAT:
	   case MPI_REAL4:
	     ierr = nc_put_vara_float (ncid, vid[nv], (size_t *) start, (size_t *) count, (const float *) bufptr); 
	     break;
	   default:
	     fprintf(stderr,"Type not recognized %d in pioc_write_darray\n",(int) iodesc->basetype);
	   }
	 }
	 break;
       case PIO_IOTYPE_NETCDF4C:
#endif
       case PIO_IOTYPE_NETCDF:
	 {
#ifdef _MPISERIAL
	   dsize = iodesc->basetype;
#else
	   mpierr = MPI_Type_size(iodesc->basetype, &dsize);
#endif
	   size_t tstart[ndims], tcount[ndims];
	   if(ios->io_rank==0){
	     for(int iorank=0;iorank<iodesc->num_aiotasks;iorank++){
	       if(iorank==0){	    
	         buflen=1;
		 for(j=0;j<ndims;j++){
		   tstart[j] =  start[j];
		   tcount[j] =  count[j];
		   buflen *= tcount[j];
		 }
	       }else{
		 mpierr = MPI_Send( &ierr, 1, MPI_INT, iorank, 0, ios->io_comm);  // handshake - tell the sending task I'm ready
		 mpierr = MPI_Recv( &buflen, 1, MPI_INT, iorank, 1, ios->io_comm, &status);
		 if(buflen>0){
		   mpierr = MPI_Recv( tstart, ndims, MPI_OFFSET, iorank, ios->num_iotasks+iorank, ios->io_comm, &status);
		   mpierr = MPI_Recv( tcount, ndims, MPI_OFFSET, iorank,2*ios->num_iotasks+iorank, ios->io_comm, &status);
		 }
	       }

	       if(buflen>0){
  	         for(int nv=0; nv<nvars; nv++){
		   if(vdesc->record >= 0){
		     tstart[0] = frame[nv];
		   }

		   //		   printf("%s %d %d %d %d %d %ld\n",__FILE__,__LINE__,iodesc->llen, buflen, iodesc->maxiobuflen, nvars, tcount[0]);


		   if(iorank>0){
		     bufptr = malloc(buflen*tsize) ;
		     mpierr = MPI_Recv( bufptr, buflen, iodesc->basetype, iorank, iorank, ios->io_comm, &status);
		   }else{
		     bufptr = (void *)((char *) IOBUF+ tsize*(nv*iodesc->llen + region->loffset));
		   }
		   //		   printf("%s %d %d %X %X\n",__FILE__,__LINE__,vid[nv],IOBUF,bufptr);
		   if(iodesc->basetype == MPI_INTEGER){
		     ierr = nc_put_vara_int (ncid, vid[nv], tstart, tcount, (const int *) bufptr); 
		   }else if(iodesc->basetype == MPI_DOUBLE || iodesc->basetype == MPI_REAL8){
		     ierr = nc_put_vara_double (ncid, vid[nv], tstart, tcount, (const double *) bufptr); 
		   }else if(iodesc->basetype == MPI_FLOAT || iodesc->basetype == MPI_REAL4){
		     ierr = nc_put_vara_float (ncid,vid[nv], tstart, tcount, (const float *) bufptr); 
		   }else{
		     fprintf(stderr,"Type not recognized %d in pioc_write_darray\n",(int) iodesc->basetype);
		   }
		   if(iorank>0){
		     free(bufptr);
		   }
		   if(ierr != PIO_NOERR){
		     for(i=0;i<fndims;i++)
		       fprintf(stderr,"vid %d dim %d start %ld count %ld\n",vid[nv],i,tstart[i],tcount[i]);
		   }
	         }
	       }
	     }
	   }else if(ios->io_rank < iodesc->num_aiotasks ){
	     buflen=1;
	     for(i=0;i<ndims;i++){
	       tstart[i] = (size_t) start[i];
	       tcount[i] = (size_t) count[i];
	       buflen*=tcount[i];
	       //	                      printf("%s %d %d %d %d\n",__FILE__,__LINE__,i,tstart[i],tcount[i]);
	     }
	     //	     printf("%s %d %d %d %d %d %d %d %d %d\n",__FILE__,__LINE__,ios->io_rank,tstart[0],tstart[1],tcount[0],tcount[1],buflen,ndims,fndims);
	     mpierr = MPI_Recv( &ierr, 1, MPI_INT, 0, 0, ios->io_comm, &status);  // task0 is ready to recieve
	     mpierr = MPI_Rsend( &buflen, 1, MPI_INT, 0, 1, ios->io_comm);
	     if(buflen>0) {
	       mpierr = MPI_Rsend( tstart, ndims, MPI_OFFSET, 0, ios->num_iotasks+ios->io_rank, ios->io_comm);
	       mpierr = MPI_Rsend( tcount, ndims, MPI_OFFSET, 0,2*ios->num_iotasks+ios->io_rank, ios->io_comm);

	       for(int nv=0; nv<nvars; nv++){
       	         bufptr = (void *)((char *) IOBUF + tsize*(nv*iodesc->llen + region->loffset));
		 //		 printf("%d %d %d\n",__LINE__,iodesc->llen,region->loffset);
		 // printf("%d %d %d %d %d %d\n",__LINE__,((int *)bufptr)[0],((int *)bufptr)[1],((int *)bufptr)[2],((int *)bufptr)[3],((int *)bufptr)[4]);
	         mpierr = MPI_Rsend( bufptr, buflen, iodesc->basetype, 0, ios->io_rank, ios->io_comm);
	       }
	     }
	   }
	   break;
	 }
	 break;
#endif
#ifdef _PNETCDF
       case PIO_IOTYPE_PNETCDF:
	 for( i=0,dsize=1;i<ndims;i++){
	   dsize*=count[i];
	 }
	 tdsize += dsize;

	 if(dsize>0){
	   //	   printf("%s %d %d %d\n",__FILE__,__LINE__,ios->io_rank,dsize);
	   startlist[rrcnt] = (PIO_Offset *) calloc(fndims, sizeof(PIO_Offset));
	   countlist[rrcnt] = (PIO_Offset *) calloc(fndims, sizeof(PIO_Offset));
	   for( i=0; i<fndims;i++){
	     startlist[rrcnt][i]=start[i];
	     countlist[rrcnt][i]=count[i];
	   }
	   rrcnt++;
	 }	 
	 if(regioncnt==iodesc->maxregions-1){
	   //printf("%s %d %d %ld %ld\n",__FILE__,__LINE__,ios->io_rank,iodesc->llen, tdsize);
	   //	   ierr = ncmpi_put_varn_all(ncid, vid, iodesc->maxregions, startlist, countlist, 
	   //			     IOBUF, iodesc->llen, iodesc->basetype);
	   
	   //printf("%s %d %ld \n",__FILE__,__LINE__,IOBUF);
	   for(int nv=0; nv<nvars; nv++){
	     if(vdesc->record >= 0){
	       for(int rc=0;rc<rrcnt;rc++){
		 startlist[rc][0] = frame[nv];
	       }
	     }
	     bufptr = (void *)((char *) IOBUF + nv*tsize*iodesc->llen);
	     ierr = ncmpi_bput_varn(ncid, vid[nv], rrcnt, startlist, countlist, 
				    bufptr, iodesc->llen, iodesc->basetype, &request);
	     pio_push_request(file,request);
	   }
	   for(i=0;i<rrcnt;i++){
	     free(startlist[i]);
	     free(countlist[i]);
	   }
	 }
	 break;
#endif
       default:
	 ierr = iotype_error(file->iotype,__FILE__,__LINE__);
       }
       if(region != NULL)
	 region = region->next;
     } //    for(regioncnt=0;regioncnt<iodesc->maxregions;regioncnt++){
   } // if(ios->ioproc)

   ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

   return ierr;
 }

int PIOc_write_darray_multi(const int ncid, const int vid[], const int ioid, const int nvars, const PIO_Offset arraylen, void *array, const int frame[], void *fillvalue)
 {
   iosystem_desc_t *ios;
   file_desc_t *file;
   io_desc_t *iodesc;
   void *iobuf;
   int vsize, rlen;
   int ierr;

   ierr = PIO_NOERR;

   file = pio_get_file_from_id(ncid);
   if(file == NULL){
     fprintf(stderr,"File handle not found %d %d\n",ncid,__LINE__);
     return PIO_EBADID;
   }
   if(! (file->mode & PIO_WRITE)){
     fprintf(stderr,"ERROR:  Attempt to write to read-only file\n");
     return PIO_EPERM;
   }
       
   iodesc = pio_get_iodesc_from_id(ioid);
   if(iodesc == NULL){
     print_trace(NULL);
     fprintf(stderr,"iodesc handle not found %d %d\n",ioid,__LINE__);
     return PIO_EBADID;
   }
   iobuf = NULL;

   pioassert(nvars>0,"nvars <= 0",__FILE__,__LINE__);

   ios = file->iosystem;
   //   rlen = iodesc->llen*nvars;
   
   rlen = iodesc->maxiobuflen*nvars;

   if(iodesc->rearranger>0){
     if(rlen>0){
       //       printf("rlen = %ld\n",rlen);
#ifdef _MPISERIAL
       vsize = iodesc->basetype;
#else
       MPI_Type_size(iodesc->basetype, &vsize);	
#endif
       iobuf = malloc( vsize * rlen);
       //     printf("%s %d 0x%.16X\n",__FILE__,__LINE__,iobuf);
/*

*/
     }
     //    printf("%s %d rlen = %d 0x%.16X 0x%.16X\n",__FILE__,__LINE__,rlen,array,iobuf); 

     //  }
     
     ierr = rearrange_comp2io(*ios, iodesc, array, iobuf, nvars);
     
   }else{
     iobuf = array;
   }
   switch(file->iotype){
   case PIO_IOTYPE_PNETCDF:
   case PIO_IOTYPE_NETCDF:
   case PIO_IOTYPE_NETCDF4P:
   case PIO_IOTYPE_NETCDF4C:
     /*     for(int ivar=0; ivar < nvars; ivar ++){
       ierr = pio_write_darray_nc(file, iodesc, vid[ivar], iobuf+vsize*ivar*(iodesc->llen), fillvalue);
     }
     */
     ierr = pio_write_darray_multi_nc(file, iodesc, nvars, vid, iobuf, frame, fillvalue);
   }
   //printf("%s %d %ld\n",__FILE__,__LINE__,iobuf);
   if(iobuf != NULL && iobuf != array){
     //     printf("%s %d 0x%.16X\n",__FILE__,__LINE__,iobuf);
     free(iobuf);
     //  printf("%s %d\n",__FILE__,__LINE__);
   }
   return ierr;

 }

 int PIOc_write_darray(const int ncid, const int vid, const int ioid, const PIO_Offset arraylen, void *array, void *fillvalue)
 {
   iosystem_desc_t *ios;
   file_desc_t *file;
   io_desc_t *iodesc;
   var_desc_t *vdesc;
   void *bufptr;
   size_t vsize, rlen;
   int ierr;
   MPI_Datatype vtype;
   wmulti_buffer *wmb;
   int tsize;
   int *tptr;
   void *bptr;
   void *fptr;
   bool recordvar;

   ierr = PIO_NOERR;

   file = pio_get_file_from_id(ncid);
   if(file == NULL){
     fprintf(stderr,"File handle not found %d %d\n",ncid,__LINE__);
     return PIO_EBADID;
   }
   if(! (file->mode & PIO_WRITE)){
     fprintf(stderr,"ERROR:  Attempt to write to read-only file\n");
     return PIO_EPERM;
   }
      
   iodesc = pio_get_iodesc_from_id(ioid);
   if(iodesc == NULL){
     fprintf(stderr,"iodesc handle not found %d %d\n",ioid,__LINE__);
     return PIO_EBADID;
   }
   ios = file->iosystem;

  vdesc = (file->varlist)+vid;
  if(vdesc == NULL)
    return PIO_EBADID;

  if(vdesc->record<0){
    recordvar=false;
  }else{
    recordvar=true;
  }


   wmb = &(file->buffer);
   if(wmb->ioid == -1){
     if(recordvar){
       wmb->ioid = ioid;
     }else{
       wmb->ioid = -(ioid);
     }
   }else{
     // separate record and non-record variables
     if(recordvar){
       while(wmb->next != NULL && wmb->ioid!=ioid){
	 if(wmb->next!=NULL)
	   wmb = wmb->next;
       }
     }else{
       while(wmb->next != NULL && wmb->ioid!= -(ioid)){
	 if(wmb->next!=NULL)
	   wmb = wmb->next;
       }
     }
   }

   if(ioid != abs(wmb->ioid) ){
     wmb->next = (wmulti_buffer *) malloc(sizeof(wmulti_buffer));
     wmb=wmb->next;
     wmb->next=NULL;
     if(recordvar){
       wmb->ioid = ioid;
     }else{
       wmb->ioid = -(ioid);
     }
     wmb->totalvars=0;
     wmb->validvars=0;
     wmb->arraylen=arraylen;
     wmb->vid=NULL;
     wmb->data=NULL;
     wmb->frame=NULL;
     wmb->fillvalue=NULL;
   }

#ifdef _MPISERIAL
    tsize = iodesc->basetype;
#else
    MPI_Type_size(iodesc->basetype, &tsize);
#endif
   // At this point wmb should be pointing to a new or existing buffer 
   // so we can add the data
   if(wmb->totalvars <= wmb->validvars){
     //     printf("%s %d %X %d %d %d\n",__FILE__,__LINE__,wmb->data,wmb->validvars,arraylen,tsize);
     bptr = (void *) realloc( wmb->data, (1+wmb->validvars)*arraylen*tsize);
     if(bptr==NULL){
       // need to flush first
       printf("%s %d %d %d %d\n",__FILE__,__LINE__,ioid,wmb->validvars,arraylen);
       PIOc_write_darray_multi(ncid, wmb->vid,  ioid, wmb->validvars, arraylen, wmb->data, wmb->frame, wmb->fillvalue);
       wmb->validvars=0;
     }else{
       wmb->data = bptr; 
       tptr = (int *) realloc( wmb->vid,sizeof(int)*( 1+wmb->validvars));
       if(tptr==NULL){
	 // bad things
	 printf("%s %d %d %d %d\n",__FILE__,__LINE__,ioid,wmb->validvars,arraylen);
       }else{
	 wmb->vid = tptr;
       }
       if(vdesc->record>=0){
	 tptr = (int *) realloc( wmb->frame,sizeof(int)*( 1+wmb->validvars));
	 if(tptr==NULL){
	   // bad things
	   printf("%s %d %d %d %d\n",__FILE__,__LINE__,ioid,wmb->validvars,arraylen);
	 }else{
	   wmb->frame = tptr;
	 }
       }
       if(iodesc->needsfill){
	 fptr = realloc( wmb->fillvalue,tsize*( 1+wmb->validvars));
	 if(fptr==NULL){
	   // bad things
	   printf("%s %d %d %d %d\n",__FILE__,__LINE__,ioid,wmb->validvars,arraylen);
	 }else{
	   wmb->fillvalue = fptr;
	 }
       }
       wmb->totalvars++;
     }
   }

   if(iodesc->needsfill){
     if(fillvalue != NULL){
       memcpy((char *) wmb->fillvalue+tsize*wmb->validvars,fillvalue, tsize); 
     }else{
       vtype = (MPI_Datatype) iodesc->basetype;
       if(vtype == MPI_INTEGER){
	 int fill = PIO_FILL_INT;
	 memcpy((char *) wmb->fillvalue+tsize*wmb->validvars, &fill, tsize); 	     
       }else if(vtype == MPI_FLOAT || vtype == MPI_REAL4){
	 float fill = PIO_FILL_FLOAT;
	 memcpy((char *) wmb->fillvalue+tsize*wmb->validvars, &fill, tsize); 
       }else if(vtype == MPI_DOUBLE || vtype == MPI_REAL8){
	 double fill = PIO_FILL_DOUBLE;
	 memcpy((char *) wmb->fillvalue+tsize*wmb->validvars, &fill, tsize); 
       }else if(vtype == MPI_CHARACTER){
	 char fill = PIO_FILL_CHAR;
	 memcpy((char *) wmb->fillvalue+tsize*wmb->validvars, &fill, tsize); 
       }else{
	 fprintf(stderr,"Type not recognized %d in pioc_write_darray\n",vtype);
       }
     }
   }
 



   wmb->arraylen = arraylen;
   wmb->vid[wmb->validvars]=vid;
   bufptr = (void *)((char *) wmb->data + arraylen*tsize*wmb->validvars);
   memcpy(bufptr, array, arraylen*tsize);

   //   printf("%s %d %d %d %d %X\n",__FILE__,__LINE__,wmb->validvars,wmb->ioid,vid,bufptr);

   if(wmb->frame!=NULL){
     wmb->frame[wmb->validvars]=vdesc->record;
   }
   wmb->validvars++;

#ifdef _PNETCDF
   PIO_Offset usage;
   if(file->iotype == PIO_IOTYPE_PNETCDF){
     // make sure we have room in the buffer ;
     if(ios->ioproc){
       ierr = ncmpi_inq_buffer_usage(ncid, &usage);
       usage += (1+wmb->validvars)*tsize*(iodesc->maxiobuflen);
     }else{
       usage=0;
     }
     MPI_Allreduce(MPI_IN_PLACE, &usage, 1,  MPI_LONG_LONG,  MPI_MAX, ios->union_comm);
     if(usage >= PIO_BUFFER_SIZE_LIMIT){	   
       //       printf("%s %d %ld %ld %ld\n",__FILE__,__LINE__,usage,PIO_BUFFER_SIZE_LIMIT,wmb->validvars*tsize*(iodesc->maxiobuflen));
       if(ios->ioproc){
	 flush_output_buffer(file);
       }
       PIOc_write_darray_multi(ncid, wmb->vid,  ioid, wmb->validvars, arraylen, wmb->data, wmb->frame, wmb->fillvalue);
       wmb->validvars=0;
     }
     
   }
#endif



   
   return ierr;

 }

int pio_read_darray_nc(file_desc_t *file, io_desc_t *iodesc, const int vid, void *IOBUF)
{
  int ierr=PIO_NOERR;
  iosystem_desc_t *ios;
  var_desc_t *vdesc;
  int ndims, fndims;
  MPI_Status status;
  int i;


  ios = file->iosystem;
  if(ios == NULL)
    return PIO_EBADID;
  
  vdesc = (file->varlist)+vid;
  
  if(vdesc == NULL)
    return PIO_EBADID;
  
  ndims = iodesc->ndims;
  ierr = PIOc_inq_varndims(file->fh, vid, &fndims);
  if(fndims==ndims) 
    vdesc->record=-1;
  
  if(ios->ioproc){
    io_region *region;
    size_t start[fndims];
    size_t count[fndims];
    size_t tmp_start[fndims];
    size_t tmp_count[fndims];
    size_t tmp_bufsize=1;
    int regioncnt;
    void *bufptr;
    int tsize;

    int rrlen=0;
    PIO_Offset *startlist[iodesc->maxregions];
    PIO_Offset *countlist[iodesc->maxregions];

    // buffer is incremented by byte and loffset is in terms of the iodessc->basetype
    // so we need to multiply by the size of the basetype
    // We can potentially allow for one iodesc to have multiple datatypes by allowing the
    // calling program to change the basetype.   
    region = iodesc->firstregion;
#ifdef _MPISERIAL
    tsize = iodesc->basetype;
#else
    MPI_Type_size(iodesc->basetype, &tsize);
#endif
    if(fndims>ndims){
      ndims++;
      if(vdesc->record<0) 
	vdesc->record=0;
    }
    for(regioncnt=0;regioncnt<iodesc->maxregions;regioncnt++){
      //            printf("%s %d %d %ld %d %d\n",__FILE__,__LINE__,regioncnt,region,fndims,ndims);
      tmp_bufsize=1;
      if(region==NULL || iodesc->llen==0){
	for(i=0;i<fndims;i++){
	  start[i] = 0;
	  count[i] = 0;
	}
	bufptr=NULL;
      }else{       
	if(regioncnt==0 || region==NULL)
	  bufptr = IOBUF;
	else
	  bufptr=(void *)((char *) IOBUF + tsize*region->loffset);
	 
	//		printf("%s %d %d %d %d\n",__FILE__,__LINE__,iodesc->llen - region->loffset, iodesc->llen, region->loffset);
	
	if(vdesc->record >= 0 && fndims>1){
	  start[0] = vdesc->record;
	  for(i=1;i<ndims;i++){
	    start[i] = region->start[i-1];
	    count[i] = region->count[i-1];
	    //	    printf("%s %d %d %ld %ld\n",__FILE__,__LINE__,i,start[i],count[i]); 
	   } 
	  if(count[1]>0)
	    count[0] = 1;
	}else{
	  // Non-time dependent array
	  for(i=0;i<ndims;i++){
	    start[i] = region->start[i];
	    count[i] = region->count[i];
	     // printf("%s %d %d %ld %ld\n",__FILE__,__LINE__,i,start[i],count[i]); 
	  }
	}
      }
       
      switch(file->iotype){
#ifdef _NETCDF
#ifdef _NETCDF4
      case PIO_IOTYPE_NETCDF4P:
	switch(iodesc->basetype){
	case MPI_DOUBLE:
	case MPI_REAL8:
	  ierr = nc_get_vara_double (file->fh, vid,start,count, bufptr); 
	  break;
	case MPI_INTEGER:
	  ierr = nc_get_vara_int (file->fh, vid, start, count,  bufptr); 
	  break;
	case MPI_FLOAT:
	case MPI_REAL4:
	  ierr = nc_get_vara_float (file->fh, vid, start,  count,  bufptr); 
	  break;
	default:
	  fprintf(stderr,"Type not recognized %d in pioc_write_darray\n",(int) iodesc->basetype);
	}	
	break;
      case PIO_IOTYPE_NETCDF4C:
#endif
      case PIO_IOTYPE_NETCDF:
	if(ios->io_rank>0){
	  tmp_bufsize=1;
	  for( i=0;i<fndims; i++){
	    tmp_start[i] = start[i];
	    tmp_count[i] = count[i];
	    tmp_bufsize *= count[i];
	  }
	  MPI_Send( tmp_count, ndims, MPI_OFFSET, 0, ios->io_rank, ios->io_comm);
	  if(tmp_bufsize > 0){
	    MPI_Send( tmp_start, ndims, MPI_OFFSET, 0, ios->io_rank, ios->io_comm);
	    //	    printf("%s %d %d\n",__FILE__,__LINE__,tmp_bufsize);
	    MPI_Recv( bufptr, tmp_bufsize, iodesc->basetype, 0, ios->io_rank, ios->io_comm, &status);
	  }
	  //	  printf("%s %d %d %d %d %d %d %d\n",__FILE__,__LINE__,regioncnt,tmp_start[1],tmp_start[2],tmp_count[1],tmp_count[2], ndims);
	}else if(ios->io_rank==0){
	  for( i=ios->num_iotasks-1; i>=0; i--){
	    if(i==0){
	      for(int k=0;k<fndims;k++)
		tmp_count[k] = count[k];
	      if(regioncnt==0 || region==NULL)
		bufptr = IOBUF;
	      else
		bufptr=(void *)((char *) IOBUF + tsize*region->loffset);
	    }else{
	      MPI_Recv(tmp_count, ndims, MPI_OFFSET, i, i, ios->io_comm, &status);
	    }
	    tmp_bufsize=1;
	    for(int j=0;j<fndims; j++){
	      tmp_bufsize *= tmp_count[j];
	    }
	    //	    printf("%s %d %d %d\n",__FILE__,__LINE__,i,tmp_bufsize);
	    if(tmp_bufsize>0){
	      if(i==0){
		for(int k=0;k<fndims;k++)
		  tmp_start[k] = start[k]; 
	      }else{
		MPI_Recv(tmp_start, ndims, MPI_OFFSET, i, i, ios->io_comm, &status);
	      }		
	      if(iodesc->basetype == MPI_DOUBLE || iodesc->basetype == MPI_REAL8){
		if(i>0)
		  bufptr = malloc(tmp_bufsize *sizeof(double));
		ierr = nc_get_vara_double (file->fh, vid, tmp_start, tmp_count, bufptr); 
	      }else if(iodesc->basetype == MPI_INTEGER){
		if(i>0)
		  bufptr = malloc(tmp_bufsize *sizeof(int));
		ierr = nc_get_vara_int (file->fh, vid, tmp_start, tmp_count,  bufptr); 	     
	      }else if(iodesc->basetype == MPI_FLOAT || iodesc->basetype == MPI_REAL4){
		if(i>0)
		  bufptr = malloc(tmp_bufsize *sizeof(float));
		ierr = nc_get_vara_float (file->fh, vid, tmp_start, tmp_count,  bufptr); 
	      }else{
		fprintf(stderr,"Type not recognized %d in pioc_write_darray\n",(int) iodesc->basetype);
	      }	
	      
	      if(ierr != PIO_NOERR){
		printf("%s %d ",__FILE__,__LINE__);
		for(int j=0;j<fndims;j++)
		  printf(" %ld %ld",tmp_start[j],tmp_count[j]);
		printf("\n");
	      }
	      
	      if(i>0){
		//    printf("%s %d %d %d\n",__FILE__,__LINE__,i,tmp_bufsize);
		MPI_Rsend(bufptr, tmp_bufsize, iodesc->basetype, i, i, ios->io_comm);
		free(bufptr);
	      }
	    }
	  }
	}
	break;
#endif
#ifdef _PNETCDF
      case PIO_IOTYPE_PNETCDF:
	{
	  tmp_bufsize=1;
	  for(int j=0;j<fndims; j++){
	    tmp_bufsize *= count[j];
	  }

	  if(tmp_bufsize>0){
             startlist[rrlen] = (PIO_Offset *) malloc(fndims * sizeof(PIO_Offset));
             countlist[rrlen] = (PIO_Offset *) malloc(fndims * sizeof(PIO_Offset));

	    for(int j=0;j<fndims; j++){
	      startlist[rrlen][j] = start[j];
	      countlist[rrlen][j] = count[j];
	      //	      printf("%s %d %d %d %d %ld %ld %ld\n",__FILE__,__LINE__,realregioncnt,iodesc->maxregions, j,start[j],count[j],tmp_bufsize);
	    }
            rrlen++;
	  }

	  if(regioncnt==iodesc->maxregions-1){
	    ierr = ncmpi_get_varn_all(file->fh, vid, rrlen, startlist, 
				      countlist, IOBUF, iodesc->llen, iodesc->basetype);
	    for(i=0;i<rrlen;i++){
	      free(startlist[i]);
	      free(countlist[i]);
	    }
	  }
	}
	break;
#endif
      default:
	ierr = iotype_error(file->iotype,__FILE__,__LINE__);
	 
      }
      if(region != NULL)
	region = region->next;
    } // for(regioncnt=0;...)
  }
  
  ierr = check_netcdf(file, ierr, __FILE__,__LINE__);

  return ierr;
}

int PIOc_read_darray(const int ncid, const int vid, const int ioid, const PIO_Offset arraylen, void *array)
{
  iosystem_desc_t *ios;
  file_desc_t *file;
  io_desc_t *iodesc;
  void *iobuf=NULL;
  size_t vsize=0, rlen=0;
  int ierr;
  MPI_Datatype vtype;
 

  file = pio_get_file_from_id(ncid);

  if(file == NULL){
    fprintf(stderr,"File handle not found %d %d\n",ncid,__LINE__);
    return PIO_EBADID;
  }
  iodesc = pio_get_iodesc_from_id(ioid);
  if(iodesc == NULL){
    fprintf(stderr,"iodesc handle not found %d %d\n",ioid,__LINE__);
    return PIO_EBADID;
  }
  ios = file->iosystem;
  if(ios->iomaster){
    rlen = iodesc->maxiobuflen;
  }else{
    rlen = iodesc->llen;
  }
  if(iodesc->rearranger > 0){
    if(ios->ioproc && rlen>0){
      vtype = (MPI_Datatype) iodesc->basetype;
      if(vtype == MPI_INTEGER){
	iobuf = malloc( rlen*sizeof(int));
      }else if(vtype == MPI_FLOAT || vtype == MPI_REAL4){
	iobuf = malloc( rlen*sizeof(float));
      }else if(vtype == MPI_DOUBLE || vtype == MPI_REAL8){
	iobuf = malloc( rlen*sizeof(double));
      }else if(vtype == MPI_CHARACTER){
	iobuf = malloc( rlen*sizeof(char));
      }else{
	fprintf(stderr,"Type not recognized %d in pioc_read_darray\n",vtype);
      }
      if(iobuf == NULL){
	fprintf(stderr,"malloc failed in pioc_read_darray %d %d\n",rlen,vtype);
	return PIO_ENOMEM;
      } 
    }
  }else{
    iobuf = array;
  }

  switch(file->iotype){
  case PIO_IOTYPE_PNETCDF:
  case PIO_IOTYPE_NETCDF:
  case PIO_IOTYPE_NETCDF4P:
  case PIO_IOTYPE_NETCDF4C:
    ierr = pio_read_darray_nc(file, iodesc, vid, iobuf);
  }
  if(iodesc->rearranger > 0){
    //	      printf("%s %d %d\n",__FILE__,__LINE__,ierr);
    ierr = rearrange_io2comp(*ios, iodesc, iobuf, array, 0, 0);
    //	      printf("%s %d %d\n",__FILE__,__LINE__,ierr);

    if(rlen>0)
      free(iobuf);
  }

  return ierr;

}

int flush_output_buffer(file_desc_t *file)
{
  var_desc_t *vardesc;
  int ierr=PIO_NOERR;
#ifdef _PNETCDF
//  if(file->nreq==0)
//    return ierr;
  int status[file->nreq];

  if(file->nreq>PIO_MAX_REQUESTS){
    fprintf(stderr,"Need to increase PIO_MAX_REQUESTS %d\n",file->nreq);
  }
     printf("%s %d \n",__FILE__,__LINE__);

  ierr = ncmpi_wait_all(file->fh,file->nreq,  file->request,status);
  for(int i=0;i<file->nreq;i++){
    file->request[i]=NC_REQ_NULL;
  }
  file->nreq = 0;
     printf("%s %d \n",__FILE__,__LINE__);

#endif
  return ierr;
}

