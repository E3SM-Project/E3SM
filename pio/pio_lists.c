#include <pio.h>
#include <pio_internal.h>
#include <string.h>
#include <stdio.h>

static io_desc_t *pio_iodesc_list=NULL;
static io_desc_t *current_iodesc=NULL;
static iosystem_desc_t *pio_iosystem_list=NULL;
static file_desc_t *pio_file_list = NULL;
static file_desc_t *current_file=NULL;

void pio_add_to_file_list(file_desc_t *file)
{
  file_desc_t *cfile;
  //  assert(file != NULL);

  file->next = NULL;
  if(pio_file_list == NULL)
    pio_file_list = file;
  else{
    for(cfile = pio_file_list; cfile->next != NULL; cfile=cfile->next);
    cfile->next = file;
  }
  current_file = file;
}


void pio_push_request(file_desc_t *file, MPI_Request request)
{
  file->request[file->nreq++] = request;
}


     
file_desc_t *pio_get_file_from_id(int ncid)
{
  file_desc_t *cfile;

  cfile = NULL;

  if(current_file != NULL && current_file->fh == ncid)
    cfile=current_file;
  for(cfile=pio_file_list; cfile != NULL; cfile=cfile->next){
    if(cfile->fh == ncid){
      current_file = cfile;
      break;
    }
  }
  return cfile;
}
  
int pio_delete_file_from_list(int ncid)
{

  file_desc_t *cfile, *pfile;

  pfile = NULL;
  for(cfile=pio_file_list; cfile != NULL; cfile=cfile->next){
    if(cfile->fh == ncid){
      if(pfile == NULL){
	pio_file_list = cfile->next;
      }else{
	pfile->next = cfile->next;
      }
      if(current_file==cfile)
	current_file = pio_file_list;
      free(cfile);
      return PIO_NOERR;
    }
    pfile = cfile;
  }
  return PIO_EBADID;
}

int pio_delete_iosystem_from_list(int piosysid)
{

  iosystem_desc_t *ciosystem, *piosystem;

  piosystem = NULL;

  //  printf(" %d iosystem_top = %ld \n",__LINE__,pio_iosystem_list);

  for(ciosystem=pio_iosystem_list; ciosystem != NULL; ciosystem=ciosystem->next){
    if(ciosystem->iosysid == piosysid){
      if(piosystem == NULL){
	pio_iosystem_list = ciosystem->next;
      }else{
	piosystem->next = ciosystem->next;
      }
      free(ciosystem);
      return PIO_NOERR;
    }
    piosystem = ciosystem;
  }
  return PIO_EBADID;
}

int pio_add_to_iosystem_list(iosystem_desc_t *ios)
{
  iosystem_desc_t *cios;
  int i=1;

  //assert(ios != NULL);
  ios->next = NULL;
  cios = pio_iosystem_list;
  if(cios==NULL)
    pio_iosystem_list = ios;
  else{
    i++;
    while(cios->next != NULL){
      cios = cios->next;
      i++;
    }
    cios->next = ios;
  }
  ios->iosysid = i << 16;
  //  ios->iosysid = i ;
  //  printf(" ios = %ld %d %ld\n",ios, ios->iosysid,ios->next);
  return ios->iosysid;
}

iosystem_desc_t *pio_get_iosystem_from_id(int iosysid)
{
  iosystem_desc_t *ciosystem;

  ciosystem = pio_iosystem_list;
  while(ciosystem != NULL){
    // printf("%d ciosystem %ld %ld %d\n",__LINE__,pio_iosystem_list, ciosystem,iosysid);
    if(ciosystem->iosysid == iosysid){
      return ciosystem;
    }
    ciosystem = ciosystem->next;
  }
  return NULL;
  
}

int pio_add_to_iodesc_list(io_desc_t *iodesc)
{
  io_desc_t *ciodesc;
  int imax=512;

  iodesc->next = NULL;
  if(pio_iodesc_list == NULL)
    pio_iodesc_list = iodesc;
  else{
    imax++;
    for(ciodesc = pio_iodesc_list; ciodesc->next != NULL; ciodesc=ciodesc->next, imax=ciodesc->ioid+1);
    ciodesc->next = iodesc;
  }
  iodesc->ioid = imax;
  current_iodesc = iodesc;
  //  printf("In add to list %d\n",iodesc->ioid);
  return iodesc->ioid;
}

     
io_desc_t *pio_get_iodesc_from_id(int ioid)
{
  io_desc_t *ciodesc;

  ciodesc = NULL;

  if(current_iodesc != NULL && current_iodesc->ioid == ioid)
    ciodesc=current_iodesc;
  for(ciodesc=pio_iodesc_list; ciodesc != NULL; ciodesc=ciodesc->next){
    if(ciodesc->ioid == ioid){
      current_iodesc = ciodesc;
      break;
    }
  }
  return ciodesc;
}
  
int pio_delete_iodesc_from_list(int ioid)
{

  io_desc_t *ciodesc, *piodesc;

  piodesc = NULL;
  for(ciodesc=pio_iodesc_list; ciodesc != NULL; ciodesc=ciodesc->next){
    if(ciodesc->ioid == ioid){
      if(piodesc == NULL){
	pio_iodesc_list = ciodesc->next;
      }else{
	piodesc->next = ciodesc->next;
      }
      if(current_iodesc==ciodesc)
	current_iodesc=pio_iodesc_list;
      free(ciodesc);
      return PIO_NOERR;
    }
    piodesc = ciodesc;
  }
  return PIO_EBADID;
}
