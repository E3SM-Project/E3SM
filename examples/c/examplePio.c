/**
 * @file 
 * @brief A simple C example for the ParallelIO Library.
 *
 * This example creates a netCDF output file with one dimension and
 * one variable. It first writes, then reads the sample file using the
 * ParallelIO library. 
 *
 * This example can be run in parallel for any number of processors up
 * to 16.
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <pio.h>
#ifdef TIMING
#include <gptl.h>
#endif

/** The length of our 1-d data array. */
static const int LEN = 16;

/** A sample data value we will use. */
static const int VAL = 42;

/** Error code in case something goes wrong. */
static const int ERR_CODE = 99;

/** Holds code and data for this example. This struct stores pointers
    to the functions used in the example, and all the data values
    needed to run the example code. */
typedef struct examplePioClass
{

    /** Pointer to initialization function. */
    struct examplePioClass* (*init)  (struct examplePioClass*);

    /** Pointer to function that creates decomposition. */
    struct examplePioClass* (*createDecomp) (struct examplePioClass*);

    /** Pointer to function that creates netCDF file. */
    struct examplePioClass* (*createFile)  (struct examplePioClass*);

    /** Pointer to function that defines netCDF variable. */
    struct examplePioClass* (*defineVar)  (struct examplePioClass*);

    /** Pointer to function that writes netCDF variable. */
    struct examplePioClass* (*writeVar) (struct examplePioClass*);

    /** Pointer to function that reads netCDF variable. */
    struct examplePioClass* (*readVar)  (struct examplePioClass*);

    /** Pointer to function that closes netCDF file. */
    struct examplePioClass* (*closeFile)  (struct examplePioClass*);

    /** Pointer to function that cleans up example memory, and library resources. */
    struct examplePioClass* (*cleanUp) (struct examplePioClass*);
    
    /** Pointer to function that handles errors. */
    struct examplePioClass* (*errorHandler) (struct examplePioClass*, const char*, const int);

    int someThing;

    /** Zero-based rank of processor. */
    int myRank;

    /** Number of processors involved in current. */
    int ntasks;

    /** Number of processors that will do IO. */
    int niotasks;

    /** Stride in the mpi rank between io tasks. */
    int stride;

    /** Number of the aggregator? */
    int numAggregator;

    /** */
    int optBase;

    /** Specifies the flavor of netCDF output format. */
    int iotype;

    /** The dimension ID. */
    int pioDimId;

    /** */
    PIO_Offset ista;

    /** */
    PIO_Offset isto;

    /** Array index per processing unit. This is the number of
     * elements of the data array that will be handled by each
     * processor. In this example there are 16 data elements. If the
     * example is run on 4 processors, then arrIdxPerPe will be 4. */
    PIO_Offset arrIdxPerPe;

    /* Length of the dimension in data. */
    int dimLen[1];

    /** The ID for the parallel I/O system. It is set by
     * PIOc_Init_Intracomm(). It references an internal structure
     * containing the general IO subsystem data and MPI structure. */
    int pioIoSystem;

    /** The ncid of the netCDF file created in this example. */
    int pioFileDesc;

    /** The ID of the netCDF varable in the example file. */
    int pioVar;

    /** The I/O description ID as passed back by PIOc_InitDecomp(). */
    int iodescNCells;

    /** A buffer for sample data. */
    int *dataBuffer;

    /** A buffer for reading data back from the file. */
    int *readBuffer;

    /** A 1-D array which holds the decomposition mapping for this
     * example. */
    PIO_Offset *compdof;

    /** The example file name. */
    char *fileName;     
    
} examplePioClass;

/** @brief Initialize libraries, create sample data. 
    
    This function is called as part of the creation of a sample data
    file for this example.

    Ths funtion initializes MPI and the ParallelIO libraries.  It sets
    up the ParallelIO library communicator. It also allocates memory
    for data used in this example, and assigns sample values to the
    data array that will be written.

    @param [in] this  Pointer to self.
    @retval examplePioClass* Pointer to self.
 */
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
    
    if (this->myRank == 0){
        printf("Running with %d MPI processes and %d PIO processes. \n",this->ntasks,this->niotasks);
    }
    
    PIOc_Init_Intracomm(MPI_COMM_WORLD, this->niotasks, this->stride, this->optBase, PIO_REARR_SUBSET, &this->pioIoSystem);
    
    /*
    ** set up some data that we will write to a netcdf file
    */
    
    this->arrIdxPerPe = LEN / this->ntasks;
    
    if (this->arrIdxPerPe < 1) {
        this->errorHandler(this, "Not enough work to distribute among pes",ERR_CODE);
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
        this->readBuffer[i] = 99;
        
        if (localVal > this->isto) {
            printf("error, should ABORT \n");
        }
        localVal++;
    }
    
    return this;
}

/** @brief Create the decomposition.

    This function is called as part of the creation of a sample data
    file for this example.

    Uses PIOc_InitDecomp() to initalize the decomposition for this
    example. The arguments are:
    - the ID of the IO system, obtained from PIOc_init_intracomm().
    - the NetCDF type of the sample data - in this case a 4-byte integer.
    - the number of dimensions (1).
    - the lengths of the dimensions.
    - the number of data elements assigned to each processor.
    - the array which provides the decomposition mapping.
    - the IO description pointer.
    - the ParallelIO rearranger (we are providing NULL here to get the
      default arranger).
    - optional array of start values for block cyclic decompositions
      (NULL means don't use block cyclical decompositions).
    - optional array of count values for block cyclic decompositions
      (NULL means don't use block cyclical decompositions).
    
    @param [in] this  Pointer to self.
    @retval examplePioClass* Pointer to self.
 */
struct examplePioClass* epc_createDecomp( struct examplePioClass* this )
{
    PIOc_InitDecomp(this->pioIoSystem, PIO_INT, 1, this->dimLen, (PIO_Offset)(this->arrIdxPerPe),
                    this->compdof, &this->iodescNCells, NULL, NULL, NULL);
    
    return this;
}

/** @brief Create the netCDF file.

    This function is called as part of the creation of a sample data
    file for this example.

    Uses the function PIOc_createfile() to create the netCDF output
    file. The format of the file is created in accordance with the
    iotype member variable, which specifies one of the following
    values: 

    - PIO_IOTYPE_PNETCDF=1 Parallel Netcdf (parallel) 
    - PIO_IOTYPE_NETCDF=2 Netcdf3 Classic format (serial) 
    - PIO_IOTYPE_NETCDF4C=3 NetCDF4 (HDF5) compressed format (serial) 
    - PIO_IOTYPE_NETCDF4P=4 NetCDF4 (HDF5) parallel

    The PIOc_createfile() function has the following parameters:

    - The IO system ID as set by PIOc_init_intracomm().
    - A pointer which will get the ncid of this file when it is created.
    - The iotype (one of the values listed above).
    - the name of the sample file.
    - the NetCDF file creating mode, PIO_CLOBBER means overwrite any
      existing file with this name.
    
    @param [in] this  Pointer to self.
    @retval examplePioClass* Pointer to self.
 */
struct examplePioClass* epc_createFile( struct examplePioClass* this )
{
    PIOc_createfile(this->pioIoSystem, &this->pioFileDesc, &this->iotype, this->fileName, PIO_CLOBBER);
    return this;
}

/** @brief Define netCDF metadata.
    
    This function is called as part of the creation of a sample data
    file for this example.

    It defines a dimension and a one-dimensional variable in the
    netCDF file using functions PIOc_def_dim() and PIOc_def_var(). It
    then calls PIOc_enddef() to end the define mode of the file.

    All of the functions take the pioFileDesc returned by
    PIOc_createfile(). This is the ncid of the netCDF file.
    
    @param [in] this  Pointer to self.
    @retval examplePioClass* Pointer to self.
 */
struct examplePioClass* epc_defineVar( struct examplePioClass* this )
{
    
    PIOc_def_dim(this->pioFileDesc, "x", (PIO_Offset)this->dimLen[0], &this->pioDimId);
    PIOc_def_var(this->pioFileDesc, "foo", PIO_INT, 1, &this->pioDimId, &this->pioVar);
    PIOc_enddef(this->pioFileDesc);

    return this;
}

/** @brief Write the sample data to the file.

    This function is called as part of the creation of a sample data
    file for this example.

    The data are written with the PIOc_write_darray() function. After
    the write is complete, ensure the file is synced for all processes
    after the write.
    
    @param [in] this  Pointer to self.
    @retval examplePioClass* Pointer to self.
 */
struct examplePioClass* epc_writeVar( struct examplePioClass* this )
{
    PIOc_write_darray(this->pioFileDesc, this->pioVar, this->iodescNCells,
                      (PIO_Offset)this->arrIdxPerPe, this->dataBuffer, NULL);
    PIOc_sync(this->pioFileDesc);
    
    return this;
}

/** @brief Read the example data from the file.

    This function is called as part of the creation of a sample data
    file for this example.

    This function reads the data that has been written to the sample
    data file. The data are read with the PIOc_read_darray() function.
    
    @param [in] this  Pointer to self.
    @retval examplePioClass* Pointer to self.
 */
struct examplePioClass* epc_readVar( struct examplePioClass* this )
{
    int i;
    
    PIOc_read_darray(this->pioFileDesc, this->pioVar, this->iodescNCells,
                     (PIO_Offset)this->arrIdxPerPe, this->readBuffer);
    
    return this;
}

/** @brief Closes the netCDF file.

    Uses the PIOc_closefile() function to close the netCDF sample file
    written by this example.

    @param [in] this Pointer to self.
    @retval examplePioClass* Pointer to self.
 */
struct examplePioClass* epc_closeFile( struct examplePioClass* this )
{
    PIOc_closefile(this->pioFileDesc);
    
    return this;
}

/** @brief Clean up allocated resources.

    This function frees the memory used in this example. It calls the
    ParallelIO library function PIOc_freedecomp() to free
    decomposition resources. Then calles PIOc_finalize() and
    MPI_finalize() to free library resources.
    
    @param [in] this  Pointer to self.
    @retval examplePioClass* Pointer to self.
 */
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

/** @brief Error handling function.

    On error, process with rank zero will print error message, the
    netCDF file will be closed with PIOc_closefile(), and MPI_Abort is
    called to end the example execution on all processes.
    
    @param [in] this  Pointer to self.
    @param [in] errMsg  an error message
    @param [in] retVal  the non-zero return value that indicated an error
    @retval examplePioClass* Pointer to self.
 */
struct examplePioClass* epc_errorHandler(struct examplePioClass* this, const char* errMsg, const int retVal)
{
   /* class(pioExampleClass), intent(inout) :: this
    character(len=*),       intent(in)    :: errMsg
    integer,                intent(in)    :: retVal*/
    
    if (retVal != PIO_NOERR){
        
        if (this->myRank == 0){
            printf("%d %s\n",retVal,errMsg);
        }
    
        PIOc_closefile(this->pioFileDesc);
        MPI_Abort(MPI_COMM_WORLD, retVal);
    }
    
    return this;
}

/** @brief Create an examplePioClass object.

    This function allocates memory for the struct that contains the
    code and data for this example. Then pointers are to the functions
    used in the example.
    
    @param [in] this  Pointer to self.
    @retval examplePioClass* Pointer to self.
 */
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
    this->errorHandler = epc_errorHandler;
    
    return this;
}

/** @brief Main execution of code.

    Executes the functions to:
    - create a new examplePioClass instance
    - initialize MPI and the ParallelIO libraries
    - create the decomposition for this example
    - create the netCDF output file
    - define the variable in the file
    - write data to the variable in the file using decomposition
    - read the data back from the file using decomposition
    - close the file
    - clean up resources

    The example can be run from the command line (on system that support it) like this:
    <pre>
    mpiexec -n 4 ./examplePio
    </pre>
    
    @param [in] argc argument count (should be zero)
    @param [in] argv argument array (should be NULL)
    @retval examplePioClass* Pointer to self.
 */
int main(int argc, const char* argv[])
{
    struct examplePioClass* pioExInst = epc_new();
    
#ifdef TIMING    
    /* Initialize the GPTL timing library. */
    int ret;
    if ((ret = GPTLinitialize ()))
      return ret;
#endif    

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
