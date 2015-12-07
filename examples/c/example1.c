/**
 * @file 
 * @brief A simple C example for the ParallelIO Library.
 *
 * This example creates a netCDF output file with one dimension and
 * one variable. It first writes, then reads the sample file using the
 * ParallelIO library. 
 *
 * This example can be run in parallel for 1, 2, 4, 8, or 16
 * processors.
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>
#include <pio.h>
#ifdef TIMING
#include <gptl.h>
#endif

#define EXAMPLE_FILENAME "example1.nc"

/* Error handling code derived from an MPI example here: 
   http://www.dartmouth.edu/~rc/classes/intro_mpi/mpi_error_functions.html */
#define MPIERR(e) do {                                                  \
	MPI_Error_string(e, err_buffer, &resultlen);			\
	printf("MPI error, line %d, file %s: %s\n", __LINE__, __FILE__, err_buffer); \
	MPI_Finalize();							\
	return 2;							\
    } while (0) 

#define ERR(e) do {				\
	MPI_Finalize();				\
	return e;				\
    } while (0) 

/* global err buffer for MPI. */
int resultlen;
char err_buffer[MPI_MAX_ERROR_STRING];

/** @brief Initialize libraries, create sample data. 
    
    This function is called as part of the creation of a sample data
    file for this example.

    Ths funtion initializes MPI and the ParallelIO libraries.  It sets
    up the ParallelIO library communicator. It also allocates memory
    for data used in this example, and assigns sample values to the
    data array that will be written.

    The ParallelIO communicator is set up with a call to
    PIOc_Init_Intracomm(). This call takes the following parameters:

    - The MPI communicator specifying the invovled processors
    (MPI_COMM_WORLD, in this case, to use all processors).
    - The number of I/O tasks. In this example there will be one I/O
    task for each process.
    - The stride (1 in this case).
    - The index of the first I/O task.
    - The iotype, specifying the flavor of netCDF to use.
    - Specify the subset rearranger.
    - A pointer that will get the ID of the ParallelIO system created
    for this call. This ID will be needed when reading or writing to
    the file using the ParallelIO library.

    @param [in] this  Pointer to self.
    @retval examplePioClass* Pointer to self.
*/

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

    The sample file created by this program is a small netCDF file. It
    has the following contents (as shown by ncdump) for a 4-processor
    run:

    <pre>
    netcdf examplePio_c {
    dimensions:
    x = 16 ;
    variables:
    int foo(x) ;
    data:

    foo = 42, 42, 42, 42, 43, 43, 43, 43, 44, 44, 44, 44, 45, 45, 45, 45 ;
    }
    </pre>
    
    @param [in] argc argument count (should be zero)
    @param [in] argv argument array (should be NULL)
    @retval examplePioClass* Pointer to self.
*/
    int main(int argc, char* argv[])
    {
	/** Set to non-zero to get output to stdout. */
	int verbose = 0;

	/** Zero-based rank of processor. */
	int my_rank;

	/** Number of processors involved in current execution. */
	int ntasks;

	/** Number of processors that will do IO. In this example we
	 * will do IO from all processors. */
	int niotasks;

	/** Stride in the mpi rank between io tasks. Always 1 in this
	 * example. */
	int stride = 1;

	/** Number of the aggregator? Always 0 in this example. */
	int numAggregator = 0;

	/** */
	int optBase = 1;

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
	int pio_io_system;

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
	char file_name[] = EXAMPLE_FILENAME;     

	/** Used for command line processing. */
	int c;

	/* Parse command line. */
	while ((c = getopt(argc, argv, "v")) != -1)
	    switch (c)
	    {
	    case 'v':
		verbose++;
		break;
	    default:
		break;
	    }

#ifdef TIMING    
	/* Initialize the GPTL timing library. */
	int ret;
	if ((ret = GPTLinitialize ()))
	    return ret;
#endif    
    
	/* Initialize MPI. */
	if ((ret = MPI_Init(&argc, &argv)))
	    MPIERR(ret);
	if ((ret = MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN)))
	    MPIERR(ret);

	/* Learn my rank and the total number of processors. */
	if ((ret = MPI_Comm_rank(MPI_COMM_WORLD, &my_rank)))
	    MPIERR(ret);
	if ((ret = MPI_Comm_size(MPI_COMM_WORLD, &ntasks)))
	    MPIERR(ret);
	/* Check that a valid number of processors was specified. */
	if (!(ntasks == 1 || ntasks == 2 || ntasks == 4 ||
	      ntasks == 8 || ntasks == 16))
	    fprintf(stderr, "Number of processors must be 1, 2, 4, 8, or 16!\n");
	if (verbose)
	    printf("%d: ParallelIO Library example1 running on %d processors.\n",
		   my_rank, ntasks);

	/* Initialize the ParallelIO library IO system. */
	iotype = PIO_IOTYPE_NETCDF;

	/* keep things simple - 1 iotask per MPI process */    
	niotasks = ntasks; 

	/* Initialize the IO system. */
	if ((ret = PIOc_Init_Intracomm(MPI_COMM_WORLD, niotasks, stride, optBase, PIO_REARR_SUBSET,
				       &pio_io_system)))
	    ERR(ret);

	/* Finalize the IO system. */
	if ((ret = PIOc_finalize(pio_io_system)))
	    ERR(ret);

	/* Finalize the MPI library. */
	MPI_Finalize();

#ifdef TIMING    
	/* Finalize the GPTL timing library. */
	if ((ret = GPTLfinalize ()))
	    return ret;
#endif    

	return 0;
    }
