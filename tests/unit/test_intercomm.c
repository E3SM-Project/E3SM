/**
 * @file Tests for PIOc_Intercomm. This tests the Init_Intercomm()
 * function, and basic asynch I/O capability.
 *
 */
#include <pio.h>
#ifdef TIMING
#include <gptl.h>
#endif

/** The number of possible output netCDF output flavors available to
 * the ParallelIO library. */
#define NUM_NETCDF_FLAVORS 4

/** The number of dimensions in the test data. */
#define NDIM 1

/** The length of our test data. */
#define DIM_LEN 4

/** The length of data on each (of 2) computational task. */
#define LOCAL_DIM_LEN 2

/** The name of the dimension in the netCDF output file. */
#define DIM_NAME "dim_test_intercomm"

/** The name of the variable in the netCDF output file. */
#define VAR_NAME "var_test_intercomm"

/** The name of the global attribute in the netCDF output file. */
#define ATT_NAME "gatt_test_intercomm"

/** Error code for when things go wrong. */
#define ERR_AWFUL 1111
#define ERR_WRONG 2222

/** Handle MPI errors. This should only be used with MPI library
 * function calls. */
#define MPIERR(e) do {                                                  \
	MPI_Error_string(e, err_buffer, &resultlen);			\
	fprintf(stderr, "MPI error, line %d, file %s: %s\n", __LINE__, __FILE__, err_buffer); \
	MPI_Finalize();							\
	return ERR_AWFUL;							\
    } while (0) 

/** Handle non-MPI errors by finalizing the MPI library and exiting
 * with an exit code. */
#define ERR(e) do {				\
        fprintf(stderr, "Error %d in %s, line %d\n", e, __FILE__, __LINE__); \
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

/** Run Tests for Init_Intercomm
 *
 * @param argc argument count
 * @param argv array of arguments
 */
int
main(int argc, char **argv)
{
    int verbose = 1;
    
    /** Zero-based rank of processor. */
    int my_rank;

    /** Number of processors involved in current execution. */
    int ntasks;

    /** Specifies the flavor of netCDF output format. */
    int iotype;

    /** Different output flavors. */
    int format[NUM_NETCDF_FLAVORS] = {PIO_IOTYPE_PNETCDF, 
				      PIO_IOTYPE_NETCDF,
				      PIO_IOTYPE_NETCDF4C,
				      PIO_IOTYPE_NETCDF4P};

    /** Names for the output files. */
    char filename[NUM_NETCDF_FLAVORS][NC_MAX_NAME + 1] = {"test_nc4_pnetcdf.nc",
							  "test_nc4_classic.nc",
							  "test_nc4_serial4.nc",
							  "test_nc4_parallel4.nc"};
	
    /** Number of processors that will do IO. In this test we
     * will do IO from all processors. */
    int niotasks;

    /** The ID for the parallel I/O system. */
    int iosysid;

    /** The ncid of the netCDF file. */
    int ncid;

    /** The ID of the netCDF varable. */
    int varid;

    /** Return code. */
    int ret;

    /** Index for loops. */
    int fmt, d, d1, i;
    
#ifdef TIMING
    /* Initialize the GPTL timing library. */
    if ((ret = GPTLinitialize ()))
	return ret;
#endif
    
    /* Initialize MPI. */
    if ((ret = MPI_Init(&argc, &argv)))
	MPIERR(ret);

    /* Learn my rank and the total number of processors. */
    if ((ret = MPI_Comm_rank(MPI_COMM_WORLD, &my_rank)))
	MPIERR(ret);
    if ((ret = MPI_Comm_size(MPI_COMM_WORLD, &ntasks)))
	MPIERR(ret);

    /* Check that a valid number of processors was specified. */
    if (!(ntasks == 1 || ntasks == 2 || ntasks == 4 ||
	  ntasks == 8 || ntasks == 16))
	fprintf(stderr, "test_intercomm Number of processors must be exactly 4!\n");
    if (verbose)
	printf("%d: test_intercomm ParallelIO Library test_intercomm running on %d processors.\n",
	       my_rank, ntasks);

    /* keep things simple - 1 iotask per MPI process */    
    niotasks = ntasks;

    /* For example, if I have 4 processors, and I want to have 2 of them be computational, */
    /* and 2 of them be IO: component count is 1  */
    /* peer_comm = MPI_COMM_WORLD */
    /* comp_comms is an array of comms of size 1 with a comm defined just over tasks (0,1) */
    /* io_comm is a comm over tasks (2,3) */

    /* Initialize the PIO IO system. This specifies how many and which
     * processors are involved in I/O. */
#define COMPONENT_COUNT 1
    MPI_Comm comp_comms;
    MPI_Comm io_comm;
    MPI_Group io_group;
    MPI_Group comp_group;

    /* Tasks 0 and 1 will be computational. Tasks 2 and 3 will be I/O
     * tasks. */

    // Get the group of processes in MPI_COMM_WORLD
    MPI_Group world_group;
    MPI_Comm_group(MPI_COMM_WORLD, &world_group);
    int comp_task;

    if (my_rank == 0 || my_rank == 1)
    {
	/* We will define comp_comm. The io_comm will get null. */
	io_comm = MPI_COMM_NULL;
	int n = 2;
	const int ranks[2] = {0, 1};

	/* Construct a group with ranks 0, 1 in world_group. */
	MPI_Group_incl(world_group, n, ranks, &comp_group);
	MPI_Comm_create_group(MPI_COMM_WORLD, comp_group, 0, &comp_comms);
	if (verbose)
	    printf("%d test_intercomm included in comp_group.\n", my_rank);

	comp_task = 1;
    }
    else
    {
	/* We will define io_comm. The comp_comms array will get nulls. */
	comp_comms = MPI_COMM_NULL;
	int n = 2;
	const int ranks[2] = {2, 3};

	/* Construct a group with ranks 2, 3 in world_group. */
	MPI_Group_incl(world_group, n, ranks, &io_group);
	MPI_Comm_create_group(MPI_COMM_WORLD, io_group, 0, &io_comm);
	if (verbose)
	    printf("%d test_intercomm included in io_group.\n", my_rank);

	comp_task = 0;
    }
    
    if ((ret = PIOc_Init_Intercomm(COMPONENT_COUNT, MPI_COMM_WORLD, &comp_comms,
				   io_comm, &iosysid)))
	ERR(ret);
    if (verbose)
	printf("%d test_intercomm init intercomm returned %d iosysid = %d\n", my_rank, ret,
	       iosysid);

    /* All the netCDF calls are only executed on the computation
     * tasks. The IO tasks have not returned from PIOc_Init_Intercomm,
     * and when the do, they should go straight to finalize. */
    if (comp_task)
    {
	for (int fmt = 0; fmt < NUM_NETCDF_FLAVORS; fmt++) 
	{
	    int ncid, varid, dimid;
	    PIO_Offset start[NDIM], count[NDIM] = {0};
	    int data[LOCAL_DIM_LEN];
	
	    /* Create a netCDF file with one dimension and one variable. */
	    if (verbose)
	    	printf("%d test_intercomm creating file %s\n", my_rank, filename[fmt]);
	    if ((ret = PIOc_createfile(iosysid, &ncid, &format[fmt], filename[fmt],
	    			       NC_CLOBBER)))
	    	ERR(ret);
	    if (verbose)
	    	printf("%d test_intercomm file created ncid = %d\n", my_rank, ncid);
	    
	    if (verbose)
	    	printf("%d test_intercomm defining dimension %s\n", my_rank, DIM_NAME);
	    if ((ret = PIOc_def_dim(ncid, DIM_NAME, DIM_LEN, &dimid)))
	    	ERR(ret);
	    if (verbose)
	    	printf("%d test_intercomm defining variable %s\n", my_rank, VAR_NAME);
	    if ((ret = PIOc_def_var(ncid, VAR_NAME, NC_INT, NDIM, &dimid, &varid)))
	    	ERR(ret);

	    /* Add a global attribute. */
	    if (verbose)
	    	printf("%d test_intercomm writing attribute %s\n", my_rank, ATT_NAME);
	    int att_data = 42;
	    /* if ((ret = PIOc_put_att_int(ncid, NC_GLOBAL, ATT_NAME, NC_INT, 1, &att_data))) */
	    /* 	ERR(ret); */

	    /* Close the file. */
	    if (verbose)
	    	printf("%d test_intercomm ending define mode ncid = %d\n", my_rank, ncid);
	    if ((ret = PIOc_enddef(ncid)))
	    	ERR(ret);

	    /* Write some data. */
	    for (int i = 0; i < LOCAL_DIM_LEN; i++)
	    	data[i] = my_rank;
	    if (verbose)
	    	printf("%d test_intercomm writing data\n", my_rank);
	    start[0] = !my_rank ? 0 : 2;
	    /* if ((ret = PIOc_put_vara_int(ncid, varid, start, count, data))) */
	    /* 	ERR(ret); */

	    /* Close the file. */
	    if (verbose)
	    	printf("%d test_intercomm closing file ncid = %d\n", my_rank, ncid);
	    if ((ret = PIOc_closefile(ncid)))
	    	ERR(ret);

	    /* Re-open the file to check it. */
	    if (verbose)
	    	printf("%d test_intercomm opening file %s\n", my_rank, filename[fmt]);
	    if ((ret = PIOc_openfile(iosysid, &ncid, &format[fmt], filename[fmt],
				     NC_NOWRITE)))
	    	ERR(ret);

	    /* Find the number of dimensions, variables, and global attributes.*/
	    int ndims, nvars, ngatts, unlimdimid;
	    if ((ret = PIOc_inq(ncid, &ndims, &nvars, &ngatts, &unlimdimid)))
	    	ERR(ret);
	    if (ndims != 1 || nvars != 1 || ngatts != 0 || unlimdimid != -1)
	    	ERR(ERR_WRONG);
	    int ndims2, nvars2, ngatts2, unlimdimid2;
	    if ((ret = PIOc_inq_ndims(ncid, &ndims2)))
	    	ERR(ret);
	    if (ndims2 != 1)
	    	ERR(ERR_WRONG);
	    if ((ret = PIOc_inq_nvars(ncid, &nvars2)))
	    	ERR(ret);
	    if (nvars2 != 1)
	    	ERR(ERR_WRONG);
	    if ((ret = PIOc_inq_natts(ncid, &ngatts2)))
	    	ERR(ret);
	    if (ngatts2 != 0)
	    	ERR(ERR_WRONG);
	    if ((ret = PIOc_inq_unlimdim(ncid, &unlimdimid2)))
	    	ERR(ret);
	    if (unlimdimid != -1)
	    	ERR(ERR_WRONG);

	    /* Check out the dimension. */
	    char dimname[NC_MAX_NAME + 1];
	    PIO_Offset dimlen;
	    if ((ret = PIOc_inq_dim(ncid, 0, dimname, &dimlen)))
	    	ERR(ret);
	    printf("%d test_intercomm dim name is %s VAR_NAME is %s\n", my_rank, dimname, VAR_NAME);
	    if (strcmp(dimname, DIM_NAME) || dimlen != DIM_LEN)
		ERR(ERR_WRONG);
	    char dimname2[NC_MAX_NAME + 1];
	    PIO_Offset dimlen2;
	    if ((ret = PIOc_inq_dimname(ncid, 0, dimname2)))
	    	ERR(ret);
	    if (strcmp(dimname2, DIM_NAME))
		ERR(ERR_WRONG);
	    if ((ret = PIOc_inq_dimlen(ncid, 0, &dimlen2)))
	    	ERR(ret);
	    if (dimlen2 != DIM_LEN)
	    	ERR(ERR_WRONG);
	    int dimid2;
	    if ((ret = PIOc_inq_dimid(ncid, DIM_NAME, &dimid2)))
	    	ERR(ret);
	    if (dimid2 != 0)
	    	ERR(ERR_WRONG);

	    /* Check out the variable. */
	    char varname[NC_MAX_NAME + 1];
	    nc_type vartype;
	    int varndims, vardimids, varnatts;
	    if ((ret = PIOc_inq_var(ncid, 0, varname, &vartype, &varndims, &vardimids, &varnatts)))
	    	ERR(ret);
	    if (strcmp(varname, VAR_NAME) || vartype != NC_INT || varndims != NDIM ||
		vardimids != 0 || varnatts != 0)
	    	ERR(ERR_WRONG);
	    char varname2[NC_MAX_NAME + 1];
	    nc_type vartype2;
	    int varndims2, vardimids2, varnatts2;
	    if ((ret = PIOc_inq_varname(ncid, 0, varname2)))
	    	ERR(ret);
	    if (strcmp(varname2, VAR_NAME))
	    	ERR(ERR_WRONG);
	    if ((ret = PIOc_inq_vartype(ncid, 0, &vartype2)))
	    	ERR(ret);
	    if (vartype2 != NC_INT)
	    	ERR(ERR_WRONG);
	    if ((ret = PIOc_inq_varndims(ncid, 0, &varndims2)))
	    	ERR(ret);
	    if (varndims2 != NDIM)
	    	ERR(ERR_WRONG);
	    if ((ret = PIOc_inq_vardimid(ncid, 0, &vardimids2)))
	    	ERR(ret);
	    if (vardimids2 != 0)
	    	ERR(ERR_WRONG);
	    if ((ret = PIOc_inq_varnatts(ncid, 0, &varnatts2)))
	    	ERR(ret);
	    if (varnatts2 != 0)
	    	ERR(ERR_WRONG);
	    int varid2;
	    if ((ret = PIOc_inq_varid(ncid, VAR_NAME, &varid2)))
	    	ERR(ret);
	    if (varid2 != 0)
	    	ERR(ERR_WRONG);
	    
	    /* Close the file. */
	    if (verbose)
	    	printf("%d test_intercomm closing file (again) ncid = %d\n", my_rank, ncid);
	    if ((ret = PIOc_closefile(ncid)))
	    	ERR(ret);

	    /* Now delete the file. */
	    if ((ret = PIOc_deletefile(iosysid, filename[fmt])))
	    	ERR(ret);
	    /* if ((ret = PIOc_openfile(iosysid, &ncid, &format[fmt], filename[fmt], */
	    /* 			     NC_NOWRITE)) != PIO_ENFILE) */
	    /* 	ERR(ERR_AWFUL); */
	    
	} /* next netcdf format flavor */
    }

    /* Free local MPI resources. */
    if (verbose)
	printf("%d test_intercomm Freeing local MPI resources...\n", my_rank);
    MPI_Group_free(&world_group);
    if (comp_task)
    {
	MPI_Group_free(&comp_group);
	MPI_Comm_free(&comp_comms);
    }
    else
    {
	MPI_Group_free(&io_group);
	MPI_Comm_free(&io_comm);
    }
    
    /* Finalize the IO system. */
    if (verbose)
	printf("%d test_intercomm Freeing PIO resources...\n", my_rank);
    if ((ret = PIOc_finalize(iosysid)))
    	ERR(ret);

    /* Finalize the MPI library. */
    MPI_Finalize();

#ifdef TIMING
    /* Finalize the GPTL timing library. */
    if ((ret = GPTLfinalize()))
	return ret;
#endif

    if (verbose)
	printf("%d test_intercomm SUCCESS!!\n", my_rank);

    
    return 0;
}
