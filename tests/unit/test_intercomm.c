/**
 * @file Tests for PIOc_Intercomm. This tests the Init_Intercomm()
 * function, and basic asynch I/O capability.
 *
 */
#include <pio.h>
#include <unistd.h>
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

/** The name of the dimension in the netCDF output file. */
#define FIRST_DIM_NAME "jojo"
#define DIM_NAME "dim_test_intercomm"

/** The name of the variable in the netCDF output file. */
#define FIRST_VAR_NAME "bill"
#define VAR_NAME "var_test_intercomm"

/** The name of the global attribute in the netCDF output file. */
#define FIRST_ATT_NAME "willy_gatt_test_intercomm"
#define ATT_NAME "gatt_test_intercomm"
#define SHORT_ATT_NAME "short_gatt_test_intercomm"
#define FLOAT_ATT_NAME "float_gatt_test_intercomm"
#define DOUBLE_ATT_NAME "double_gatt_test_intercomm"

/** The value of the global attribute in the netCDF output file. */
#define ATT_VALUE 42

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

/* Check the file for correctness. */
int
check_file(int iosysid, int format, char *filename, int my_rank, int verbose)
{
    int ncid;
    int ret;
    int ndims, nvars, ngatts, unlimdimid;
    int ndims2, nvars2, ngatts2, unlimdimid2;
    int dimid2;
    char dimname[NC_MAX_NAME + 1];
    PIO_Offset dimlen;
    char dimname2[NC_MAX_NAME + 1];
    PIO_Offset dimlen2;
    char varname[NC_MAX_NAME + 1];
    nc_type vartype;
    int varndims, vardimids, varnatts;
    char varname2[NC_MAX_NAME + 1];
    nc_type vartype2;
    int varndims2, vardimids2, varnatts2;
    int varid2;
    int att_data;
    short short_att_data;
    float float_att_data;
    double double_att_data;
        
    /* Re-open the file to check it. */
    if (verbose)
	printf("%d test_intercomm opening file %s format %d\n", my_rank, filename, format);
    if ((ret = PIOc_openfile(iosysid, &ncid, &format, filename,
			     NC_NOWRITE)))
	ERR(ret);
    
    /* Try to read the data. */
    PIO_Offset start[NDIM] = {0}, count[NDIM] = {DIM_LEN};    
    int data_in[DIM_LEN];
    if ((ret = PIOc_get_vars_tc(ncid, 0, start, count, NULL, NC_INT, data_in)))
	ERR(ret);
    for (int i = 0; i < DIM_LEN; i++)
    {
	if (verbose)
	    printf("%d test_intercomm read data_in[%d] = %d\n", my_rank, i, data_in[i]);
	if (data_in[i] != i)
	    ERR(ERR_AWFUL);
    }

    /* Find the number of dimensions, variables, and global attributes.*/
    if ((ret = PIOc_inq(ncid, &ndims, &nvars, &ngatts, &unlimdimid)))
    	ERR(ret);
    if (ndims != 1 || nvars != 1 || ngatts != 4 || unlimdimid != -1)
    	ERR(ERR_WRONG);

    /* This should return PIO_NOERR. */
    if ((ret = PIOc_inq(ncid, NULL, NULL, NULL, NULL)))
    	ERR(ret);

    /* Check the other functions that get these values. */
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
    if (ngatts2 != 4)
    	ERR(ERR_WRONG);
    if ((ret = PIOc_inq_unlimdim(ncid, &unlimdimid2)))
    	ERR(ret);
    if (unlimdimid != -1)
    	ERR(ERR_WRONG);
    /* Should succeed, do nothing. */   
    if ((ret = PIOc_inq_unlimdim(ncid, NULL))) 
    	ERR(ret);

    /* Check out the dimension. */
    if ((ret = PIOc_inq_dim(ncid, 0, dimname, &dimlen)))
    	ERR(ret);
    if (strcmp(dimname, DIM_NAME) || dimlen != DIM_LEN)
    	ERR(ERR_WRONG);

    /* Check the other functions that get these values. */
    if ((ret = PIOc_inq_dimname(ncid, 0, dimname2)))
    	ERR(ret);
    if (strcmp(dimname2, DIM_NAME))
    	ERR(ERR_WRONG);
    if ((ret = PIOc_inq_dimlen(ncid, 0, &dimlen2)))
    	ERR(ret);
    if (dimlen2 != DIM_LEN)
    	ERR(ERR_WRONG);
    if ((ret = PIOc_inq_dimid(ncid, DIM_NAME, &dimid2)))
    	ERR(ret);
    if (dimid2 != 0)
    	ERR(ERR_WRONG);

    /* Check out the variable. */
    if ((ret = PIOc_inq_var(ncid, 0, varname, &vartype, &varndims, &vardimids, &varnatts)))
    	ERR(ret);
    if (strcmp(varname, VAR_NAME) || vartype != NC_INT || varndims != NDIM ||
    	vardimids != 0 || varnatts != 0)
    	ERR(ERR_WRONG);

    /* Check the other functions that get these values. */
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
    if ((ret = PIOc_inq_varid(ncid, VAR_NAME, &varid2)))
    	ERR(ret);
    if (varid2 != 0)
    	ERR(ERR_WRONG);

    /* Check out the global attributes. */
    nc_type atttype;
    PIO_Offset attlen;
    char myattname[NC_MAX_NAME + 1];
    int myid;
    if ((ret = PIOc_inq_att(ncid, NC_GLOBAL, ATT_NAME, &atttype, &attlen)))
    	ERR(ret);
    if (atttype != NC_INT || attlen != 1)
    	ERR(ERR_WRONG);
    if ((ret = PIOc_inq_attlen(ncid, NC_GLOBAL, ATT_NAME, &attlen)))
    	ERR(ret);
    if (attlen != 1)
    	ERR(ERR_WRONG);
    if ((ret = PIOc_inq_attname(ncid, NC_GLOBAL, 0, myattname)))
    	ERR(ret);
    if (strcmp(ATT_NAME, myattname))
    	ERR(ERR_WRONG);
    if ((ret = PIOc_inq_attid(ncid, NC_GLOBAL, ATT_NAME, &myid)))
    	ERR(ret);
    if (myid != 0)
    	ERR(ERR_WRONG);
    if ((ret = PIOc_get_att_int(ncid, NC_GLOBAL, ATT_NAME, &att_data)))
    	ERR(ret);
    if (verbose)
    	printf("%d test_intercomm att_data = %d\n", my_rank, att_data);
    if (att_data != ATT_VALUE)
    	ERR(ERR_WRONG);
    if ((ret = PIOc_inq_att(ncid, NC_GLOBAL, SHORT_ATT_NAME, &atttype, &attlen)))
    	ERR(ret);
    if (atttype != NC_SHORT || attlen != 1)
    	ERR(ERR_WRONG);
    if ((ret = PIOc_get_att_short(ncid, NC_GLOBAL, SHORT_ATT_NAME, &short_att_data)))
    	ERR(ret);
    if (short_att_data != ATT_VALUE)
    	ERR(ERR_WRONG);
    if ((ret = PIOc_get_att_float(ncid, NC_GLOBAL, FLOAT_ATT_NAME, &float_att_data)))
    	ERR(ret);
    if (float_att_data != ATT_VALUE)
    	ERR(ERR_WRONG);
    if ((ret = PIOc_get_att_double(ncid, NC_GLOBAL, DOUBLE_ATT_NAME, &double_att_data)))
    	ERR(ret);
    if (double_att_data != ATT_VALUE)
    	ERR(ERR_WRONG);
    
	    
    /* Close the file. */
    if (verbose)
	printf("%d test_intercomm closing file (again) ncid = %d\n", my_rank, ncid);
    if ((ret = PIOc_closefile(ncid)))
	ERR(ret);

    return 0;
}

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

    /** Different output flavors. */
    int format[NUM_NETCDF_FLAVORS] = {PIO_IOTYPE_PNETCDF, 
				      PIO_IOTYPE_NETCDF,
				      PIO_IOTYPE_NETCDF4C,
				      PIO_IOTYPE_NETCDF4P};

    /** Names for the output files. */
    char filename[NUM_NETCDF_FLAVORS][NC_MAX_NAME + 1] = {"test_intercomm_pnetcdf.nc",
							  "test_intercomm_classic.nc",
							  "test_intercomm_serial4.nc",
							  "test_intercomm_parallel4.nc"};
	
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
    if (ntasks != 4)
    {
	fprintf(stderr, "test_intercomm Number of processors must be exactly 4!\n");
	ERR(ERR_AWFUL);
    }
    if (verbose)
	printf("%d: test_intercomm ParallelIO Library test_intercomm running on %d processors.\n",
	       my_rank, ntasks);

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

    /* Tasks 0 and 1 will be computational. Tasks 2 and 3 will be I/O
     * tasks. */
    int comp_task;

    int color = my_rank < 2 ? 0 : 1; // Determine color based on row

    // Split the communicator based on the color and use the
    // original rank for ordering
    MPI_Comm row_comm;
    MPI_Comm_split(MPI_COMM_WORLD, color, my_rank, &row_comm);

    int row_rank, row_size;
    MPI_Comm_rank(row_comm, &row_rank);
    MPI_Comm_size(row_comm, &row_size);

    printf("WORLD RANK: %d \t ROW RANK/SIZE: %d/%d\n",
	   my_rank, row_rank, row_size);
    if (my_rank == 0 || my_rank == 1)
    {
	/* We will define comp_comm. The io_comm will get null. */
	io_comm = MPI_COMM_NULL;
	comp_comms = row_comm;
	comp_task = 1;
	if (verbose)
	    printf("%d added to the comp_comm\n", my_rank);
    }
    else
    {
	/* We will define io_comm. The comp_comms array will get nulls. */
	comp_comms = MPI_COMM_NULL;
	io_comm = row_comm;
	comp_task = 0;
	if (verbose)
	    printf("%d added to the io_comm\n", my_rank);
    }

    /* Turn on logging. */
    if ((ret = PIOc_set_log_level(3)))
	ERR(ret);

    /* Initialize the async setup. */
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
	    int data[DIM_LEN];
	
	    /* Create a netCDF file with one dimension and one variable. */
	    if (verbose)
	    	printf("%d test_intercomm creating file %s\n", my_rank, filename[fmt]);
	    if ((ret = PIOc_createfile(iosysid, &ncid, &format[fmt], filename[fmt],
	    			       NC_CLOBBER)))
	    	ERR(ret);
	    if (verbose)
	    	printf("%d test_intercomm file created ncid = %d\n", my_rank, ncid);

	    /* /\* End define mode, then re-enter it. *\/ */
	    if ((ret = PIOc_enddef(ncid)))
	    	ERR(ret);
	    if (verbose)
	    	printf("%d test_intercomm calling redef\n", my_rank);
	    if ((ret = PIOc_redef(ncid)))
	    	ERR(ret);

	    /* Test the inq_format function. */
	    int myformat;
	    if ((ret = PIOc_inq_format(ncid, &myformat))) 
		ERR(ret); 
	    if ((format[fmt] == PIO_IOTYPE_PNETCDF || format[fmt] == PIO_IOTYPE_NETCDF) &&
		myformat != 1)
		ERR(ERR_AWFUL);
	    else if ((format[fmt] == PIO_IOTYPE_NETCDF4C || format[fmt] == PIO_IOTYPE_NETCDF4P) &&
		     myformat != 3)
		ERR(ERR_AWFUL);

	    /* Test the inq_type function for atomic types. */
	    char type_name[NC_MAX_NAME + 1];
	    PIO_Offset type_size;
	    #define NUM_TYPES 11
	    nc_type xtype[NUM_TYPES] = {NC_CHAR, NC_BYTE, NC_SHORT, NC_INT, NC_FLOAT, NC_DOUBLE,
					NC_UBYTE, NC_USHORT, NC_UINT, NC_INT64, NC_UINT64};
	    int type_len[NUM_TYPES] = {1, 1, 2, 4, 4, 8, 1, 2, 4, 8, 8};
	    int max_type = format[fmt] == PIO_IOTYPE_NETCDF ? NC_DOUBLE : NC_UINT64;
	    for (int i = 0; i < max_type; i++)
	    {
		if ((ret = PIOc_inq_type(ncid, xtype[i], type_name, &type_size)))
		    ERR(ret);
		if (type_size != type_len[i])
		    ERR(ERR_AWFUL);
	    }
	    
	    /* Define a dimension. */
	    char dimname2[NC_MAX_NAME + 1];
	    if (verbose)
	    	printf("%d test_intercomm defining dimension %s\n", my_rank, DIM_NAME);
	    if ((ret = PIOc_def_dim(ncid, FIRST_DIM_NAME, DIM_LEN, &dimid)))
	    	ERR(ret);
	    if ((ret = PIOc_inq_dimname(ncid, 0, dimname2)))
		ERR(ret);
	    if (strcmp(dimname2, FIRST_DIM_NAME))
		ERR(ERR_WRONG);
	    if ((ret = PIOc_rename_dim(ncid, 0, DIM_NAME)))
	    	ERR(ret);

	    /* Define a 1-D variable. */
	    char varname2[NC_MAX_NAME + 1];
	    if (verbose)
	    	printf("%d test_intercomm defining variable %s\n", my_rank, VAR_NAME);
	    if ((ret = PIOc_def_var(ncid, FIRST_VAR_NAME, NC_INT, NDIM, &dimid, &varid)))
	    	ERR(ret);
	    if ((ret = PIOc_inq_varname(ncid, 0, varname2)))
		ERR(ret);
	    if (strcmp(varname2, FIRST_VAR_NAME))
		ERR(ERR_WRONG);
	    if ((ret = PIOc_rename_var(ncid, 0, VAR_NAME)))
	    	ERR(ret);

	    /* Add a global attribute. */
	    if (verbose)
	    	printf("%d test_intercomm writing attributes %s\n", my_rank, ATT_NAME);
	    int att_data = ATT_VALUE;
	    short short_att_data = ATT_VALUE;
	    float float_att_data = ATT_VALUE;
	    double double_att_data = ATT_VALUE;
	    char attname2[NC_MAX_NAME + 1];
	    /* Write an att and rename it. */
	    if ((ret = PIOc_put_att_int(ncid, NC_GLOBAL, FIRST_ATT_NAME, NC_INT, 1, &att_data)))
		ERR(ret);
	    if ((ret = PIOc_inq_attname(ncid, NC_GLOBAL, 0, attname2)))
		ERR(ret);
	    if (strcmp(attname2, FIRST_ATT_NAME))
		ERR(ERR_WRONG);
	    if ((ret = PIOc_rename_att(ncid, NC_GLOBAL, FIRST_ATT_NAME, ATT_NAME)))
		ERR(ret);

	    /* Write an att and delete it. */
	    nc_type myatttype;
	    if ((ret = PIOc_put_att_int(ncid, NC_GLOBAL, FIRST_ATT_NAME, NC_INT, 1, &att_data)))
		ERR(ret);
	    if ((ret = PIOc_del_att(ncid, NC_GLOBAL, FIRST_ATT_NAME)))
		ERR(ret);
	    /* if ((ret = PIOc_inq_att(ncid, NC_GLOBAL, FIRST_ATT_NAME, NULL, NULL)) != PIO_ENOTATT) */
	    /* { */
	    /* 	printf("ret = %d\n", ret); */
	    /* 	ERR(ERR_AWFUL); */
	    /* } */

	    /* Write some atts of different types. */
	    if ((ret = PIOc_put_att_short(ncid, NC_GLOBAL, SHORT_ATT_NAME, NC_SHORT, 1, &short_att_data)))
	    	ERR(ret);
	    if ((ret = PIOc_put_att_float(ncid, NC_GLOBAL, FLOAT_ATT_NAME, NC_FLOAT, 1, &float_att_data)))
	    	ERR(ret);
	    if ((ret = PIOc_put_att_double(ncid, NC_GLOBAL, DOUBLE_ATT_NAME, NC_DOUBLE, 1, &double_att_data)))
	    	ERR(ret);

	    /* End define mode. */
	    if (verbose)
	    	printf("%d test_intercomm ending define mode ncid = %d\n", my_rank, ncid);
	    if ((ret = PIOc_enddef(ncid)))
	    	ERR(ret);

	    /* Write some data. For the PIOc_put/get functions, all
	     * data must be on compmaster before the function is
	     * called. Only compmaster's arguments are passed to the
	     * async msg handler. All other computation tasks are
	     * ignored. */
	    for (int i = 0; i < DIM_LEN; i++)
	    	data[i] = i;
	    if (verbose)
	    	printf("%d test_intercomm writing data\n", my_rank);
	    if (verbose)
	    	printf("%d test_intercomm writing data\n", my_rank);
	    start[0] = 0;
	    count[0] = DIM_LEN;
	    if ((ret = PIOc_put_vars_tc(ncid, varid, start, count, NULL, NC_INT, data)))
	    	ERR(ret);

	    /* Close the file. */
	    if (verbose)
	    	printf("%d test_intercomm closing file ncid = %d\n", my_rank, ncid);
	    if ((ret = PIOc_closefile(ncid)))
	    	ERR(ret);

	    /* Check the file for correctness. */
	    if ((ret = check_file(iosysid, format[fmt], filename[fmt], my_rank, verbose)))
	    	ERR(ret);

	    /* Now delete the file. */
	    /* if ((ret = PIOc_deletefile(iosysid, filename[fmt]))) */
	    /* 	ERR(ret); */
	    /* if ((ret = PIOc_openfile(iosysid, &ncid, &format[fmt], filename[fmt], */
	    /* 			     NC_NOWRITE)) != PIO_ENFILE) */
	    /* 	ERR(ERR_AWFUL); */
	    
	} /* next netcdf format flavor */
    }

    /* Free local MPI resources. */
    if (verbose)
	printf("%d test_intercomm Freeing local MPI resources...\n", my_rank);
    if (comp_task)
    {
	MPI_Comm_free(&comp_comms);
    }
    else
    {
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
