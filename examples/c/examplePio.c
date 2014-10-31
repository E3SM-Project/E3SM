#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <pio.h>

static const int LEN = 16;
static const int VAL = 42;
static const int ERR_CODE = 99;

typedef struct examplePioClass
{
    
    struct examplePioClass* (*init)  (struct examplePioClass*);
    struct examplePioClass* (*createDecomp) (struct examplePioClass*);
    struct examplePioClass* (*createFile)  (struct examplePioClass*);
    struct examplePioClass* (*defineVar)  (struct examplePioClass*);
    struct examplePioClass* (*writeVar) (struct examplePioClass*);
    struct examplePioClass* (*readVar)  (struct examplePioClass*);
    struct examplePioClass* (*closeFile)  (struct examplePioClass*);
    struct examplePioClass* (*cleanUp) (struct examplePioClass*);

    int someThing;
    
    int myRank;
    int ntasks;
    
    int niotasks;
    int stride; /* stride in the mpi rank between io tasks */
    int numAggregator;
    int optBase;
    int iotype;
    int pioDimId;
    PIO_Offset ista;
    PIO_Offset isto;
    PIO_Offset arrIdxPerPe;
    int dimLen[1];
    
    int pioIoSystem;
    int pioFileDesc;
    int pioVar;
    int iodescNCells;
    
    int *dataBuffer;
    int *readBuffer;
    PIO_Offset *compdof;
    
    char *fileName;
    
} examplePioClass;

struct examplePioClass* epc_init( struct examplePioClass* this )
{
    int ierr;
    int argc;
    char *argv;
    int i, localVal;

    /*
    ** initialize MPI
    */
    
    ierr = MPI_Init(NULL, NULL);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &this->myRank);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &this->ntasks);
    
    printf(" examplePioClass::init %d\n",this->myRank);

    /*
    ** set up PIO for rest of example
    */
        
    this->stride        = 1;
    this->numAggregator = 0;
    this->optBase       = 1;
    this->iotype        = PIO_IOTYPE_NETCDF;
    this->fileName      = "examplePio_c.nc";
    this->dimLen[0]     = LEN;
    
    this->niotasks = this->ntasks; /* keep things simple - 1 iotask per MPI process */
    
    PIOc_Init_Intracomm(MPI_COMM_WORLD, this->niotasks, this->stride, this->optBase, PIO_REARR_NONE, &this->pioIoSystem);
    
    /*
    ** set up some data that we will write to a netcdf file
    */
    
    this->arrIdxPerPe = LEN / this->ntasks;
    
    if (this->arrIdxPerPe < 1) {
        printf("not enough work to do \n");
        /*call this%errorHandle("Not enough work to distribute among pes", ERR_CODE)*/
    }
        
    this->ista = this->myRank * this->arrIdxPerPe;
    this->isto = this->ista + (this->arrIdxPerPe - 1);
    
    this->dataBuffer = (int *)malloc(this->arrIdxPerPe * sizeof (int));
    this->readBuffer = (int *)malloc(this->arrIdxPerPe * sizeof (int));
    this->compdof = (PIO_Offset *)malloc(this->arrIdxPerPe * sizeof(PIO_Offset));

    /*
    ** assign values to various arrays
    */
    
    localVal = this->ista;
    for (i = 0; i < this->arrIdxPerPe; i++ ){
        
        this->dataBuffer[i] = this->myRank + VAL;
        this->compdof[i] = localVal;
        this->readBuffer[i] = 0;
        
        if (localVal > this->isto) {
            printf("error, should ABORT \n");
        }
        localVal++;
    }
    
    return this;
}

struct examplePioClass* epc_createDecomp( struct examplePioClass* this )
{
    PIOc_InitDecomp(this->pioIoSystem, PIO_INT, 1, this->dimLen, (PIO_Offset)(this->arrIdxPerPe),
                    this->compdof, &this->iodescNCells, NULL, &this->ista, &this->arrIdxPerPe);
    return this;
}

struct examplePioClass* epc_createFile( struct examplePioClass* this )
{
    PIOc_createfile(this->pioIoSystem, &this->pioFileDesc, &this->iotype, this->fileName, PIO_CLOBBER);
    return this;
}

struct examplePioClass* epc_defineVar( struct examplePioClass* this )
{
    
    PIOc_def_dim(this->pioFileDesc, "x", (PIO_Offset)this->dimLen[0], &this->pioDimId);
    
    PIOc_def_var(this->pioFileDesc, "foo", PIO_INT, 1, &this->pioDimId, &this->pioVar);
    
    PIOc_enddef(this->pioFileDesc);

    return this;
}

struct examplePioClass* epc_writeVar( struct examplePioClass* this )
{
    PIOc_write_darray(this->pioFileDesc, this->pioVar, this->iodescNCells,
                      (PIO_Offset)this->arrIdxPerPe, this->dataBuffer, NULL);
    PIOc_sync(this->pioFileDesc);
    return this;
}

struct examplePioClass* epc_readVar( struct examplePioClass* this )
{
    int i;
    PIOc_read_darray(this->pioFileDesc, this->pioVar, this->iodescNCells,
                     (PIO_Offset)this->arrIdxPerPe, this->readBuffer);
    
    for (i = 0; i < this->arrIdxPerPe; i++ ){
        printf("after read %d %d\n", this->myRank,this->readBuffer[i]);
    }
    
    return this;
}

struct examplePioClass* epc_closeFile( struct examplePioClass* this )
{
    PIOc_closefile(this->pioFileDesc);
    return this;
}

struct examplePioClass* epc_cleanUp( struct examplePioClass* this )
{
    int ierr;
    
    free(this->dataBuffer);
    free(this->readBuffer);
    free(this->compdof);
    
    ierr = PIOc_freedecomp(this->pioIoSystem, this->iodescNCells);
    ierr = PIOc_finalize(this->pioIoSystem);
    ierr = MPI_Finalize();
    
    return this;
}

struct examplePioClass* epc_new()
{
    struct examplePioClass* this = malloc((sizeof(struct examplePioClass)));
    
    /* assign function pointers to impls */
    
    this->init = epc_init;
    this->createDecomp = epc_createDecomp;
    this->createFile = epc_createFile;
    this->defineVar = epc_defineVar;
    this->writeVar = epc_writeVar;
    this->readVar = epc_readVar;
    this->closeFile = epc_closeFile;
    this->cleanUp = epc_cleanUp;
    
    return this;
}

int main(int argc, const char* argv[])
{
    struct examplePioClass* pioExInst = epc_new();
    
    pioExInst->init(pioExInst);
    pioExInst->createDecomp(pioExInst);
    pioExInst->createFile(pioExInst);
    pioExInst->defineVar(pioExInst);
    pioExInst->writeVar(pioExInst);
    pioExInst->readVar(pioExInst);
    pioExInst->closeFile(pioExInst);
    pioExInst->cleanUp(pioExInst);
    
    free(pioExInst);
    return 0;
}