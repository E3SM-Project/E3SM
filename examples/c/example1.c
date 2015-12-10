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

/** This is the fiename used for the netCDF file created by this
 * example. */
#define FILENAME "example1.nc"

/** The number of dimensions in the example data. In this simple
    example, we are using one-dimensional data. */
#define NDIM 1

/** The length of our sample data. There will be a total of 16
 * integers in our data, and responsibilty for writing and reading
 * them will be spread between all the processors used to run this
 * example. */
#define DIM_LEN 16

/** Handle MPI errors. This should only be used with MPI library
 * function calls. */
#define MPIERR(e) do {                                                  \
	MPI_Error_string(e, err_buffer, &resultlen);			\
	printf("MPI error, line %d, file %s: %s\n", __LINE__, __FILE__, err_buffer); \
	MPI_Finalize();							\
	return 2;							\
    } while (0) 

/** Handle non-MPI errors by finalizing the MPI library and exiting
 * with an exit code. */
#define ERR(e) do {				\
	MPI_Finalize();				\
	return e;				\
    } while (0) 

/** Global err buffer for MPI. When there is an MPI error, this buffer
 * is used to store the error message that is associated with the MPI
 * error. */
char err_buffer[MPI_MAX_ERROR_STRING];

/** This is the length of the most recent MPI error message, stored
 * int the global error string. */
int resultlen;

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
	int ioproc_stride = 1;

	/** Number of the aggregator? Always 0 in this example. */
	int numAggregator = 0;

	/** Zero based rank of first processor to be used for I/O. */
	int ioproc_start = 0;

	/** Specifies the flavor of netCDF output format. */
	int iotype;

	/** The dimension ID. */
	int dimid;

	/** Array index per processing unit. This is the number of
	 * elements of the data array that will be handled by each
	 * processor. In this example there are 16 data elements. If the
	 * example is run on 4 processors, then arrIdxPerPe will be 4. */
	PIO_Offset elements_per_pe;

	/* Length of the dimensions in the data. This simple example
	 * uses one-dimensional data. The lenght along that dimension
	 * is DIM_LEN (16). */
	int dim_len[1] = {DIM_LEN};

	/** The ID for the parallel I/O system. It is set by
	 * PIOc_Init_Intracomm(). It references an internal structure
	 * containing the general IO subsystem data and MPI
	 * structure. It is passed to PIOc_finalize() to free
	 * associated resources, after all I/O, but before
	 * MPI_Finalize is called. */
	int iosysid;

	/** The ncid of the netCDF file created in this example. */
	int ncid;

	/** The ID of the netCDF varable in the example file. */
	int varid;

	/** The I/O description ID as passed back by PIOc_InitDecomp()
	 * and freed in PIOc_freedecomp(). */
	int ioid;

	/** A buffer for sample data.  The size of this array will
	 * vary depending on how many processors are involved in the
	 * execution of the example code. It's length will be the same
	 * as elements_per_pe.*/
	int *buffer;

	/** A buffer for reading data back from the file. The size of
	 * this array will vary depending on how many processors are
	 * involved in the execution of the example code. It's length
	 * will be the same as elements_per_pe.*/
	int *read_buffer;

	/** A 1-D array which holds the decomposition mapping for this
	 * example. The size of this array will vary depending on how
	 * many processors are involved in the execution of the
	 * example code. It's length will be the same as
	 * elements_per_pe. */
	PIO_Offset *compdof;

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

	/* keep things simple - 1 iotask per MPI process */    
	niotasks = ntasks; 

	/* Initialize the PIO IO system. This specifies how
	 * many and which processors are involved in I/O. */
	if ((ret = PIOc_Init_Intracomm(MPI_COMM_WORLD, niotasks, ioproc_stride,
				       ioproc_start, PIO_REARR_SUBSET, &iosysid)))
	    ERR(ret);

	/* Describe the decomposition. This is a 1-based array, so add 1! */
	elements_per_pe = DIM_LEN / ntasks;
	if (!(compdof = malloc(elements_per_pe * sizeof(PIO_Offset))))
	    return PIO_ENOMEM;
	for (int i = 0; i < elements_per_pe; i++)
	    compdof[i] = my_rank * elements_per_pe + i + 1;
	
	/* Create the PIO decomposition for this example. */
	if (verbose)
	    printf("rank: %d Creating decomposition...\n", my_rank);
	if ((ret = PIOc_InitDecomp(iosysid, PIO_INT, NDIM, dim_len, (PIO_Offset)elements_per_pe,
				   compdof, &ioid, NULL, NULL, NULL)))
	    ERR(ret);
	free(compdof);
	
	/* Set the format of the output file. */
	iotype = PIO_IOTYPE_NETCDF;
	
	/* Create the netCDF output file. */
	if (verbose)
	    printf("rank: %d Creating sample file %s with iotype %d...\n",
		   my_rank, FILENAME, iotype);
	if ((ret = PIOc_createfile(iosysid, &ncid, &iotype, FILENAME, PIO_CLOBBER)))
	    ERR(ret);
	
	/* Define netCDF dimension and variable. */
	if (verbose)
	    printf("rank: %d Defining netCDF metadata...\n", my_rank);
	if ((ret = PIOc_def_dim(ncid, "x", (PIO_Offset)dim_len[0], &dimid)))
	    ERR(ret);
	if ((ret = PIOc_def_var(ncid, "foo", PIO_INT, NDIM, &dimid, &varid)))
	    ERR(ret);
	if ((ret = PIOc_enddef(ncid)))
	    ERR(ret);
	
	/* /\* Prepare sample data. *\/ */
	/* if (!(buffer = malloc(elements_per_pe * sizeof(int)))) */
	/*     return PIO_ENOMEM; */
	/* for (int i = 0; i < elements_per_pe; i++) */
	/*     buffer[i] = 42 + my_rank; */

	/* /\* Write data to the file. *\/ */
	/* if (verbose) */
	/*     printf("rank: %d Writing sample data...\n", my_rank); */
	/* if ((ret = PIOc_write_darray(ncid, varid, ioid, (PIO_Offset)elements_per_pe, */
	/* 			     buffer, NULL))) */
	/*     ERR(ret); */
	/* if ((ret = PIOc_sync(ncid))) */
	/*     ERR(ret); */

	/* /\* Free buffer space used in this example. *\/ */
	/* free(buffer); */
	
	/* Close the netCDF file. */
	if (verbose)
	    printf("rank: %d Closing the sample data file...\n", my_rank);
	if ((ret = PIOc_closefile(ncid)))
	    ERR(ret);
	
	/* Free the PIO decomposition. */
	if (verbose)
	    printf("rank: %d Freeing PIO decomposition...\n", my_rank);
	if ((ret = PIOc_freedecomp(iosysid, ioid)))
	    ERR(ret);
	     
	
	/* Finalize the IO system. */
	if (verbose)
	    printf("rank: %d Freeing PIO resources...\n", my_rank);
	if ((ret = PIOc_finalize(iosysid)))
	    ERR(ret);

	/* Finalize the MPI library. */
	MPI_Finalize();

#ifdef TIMING    
	/* Finalize the GPTL timing library. */
	if ((ret = GPTLfinalize ()))
	    return ret;
#endif    

	if (verbose)
	    printf("rank: %d SUCCESS!\n", my_rank);
	return 0;
    }
