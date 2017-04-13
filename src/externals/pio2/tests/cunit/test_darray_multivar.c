/*
 * Tests for PIO distributed arrays.
 *
 * Ed Hartnett, 2/16/17
 */
#include <pio.h>
#include <pio_internal.h>
#include <pio_tests.h>

/* The number of tasks this test should run on. */
#define TARGET_NTASKS 4

/* The minimum number of tasks this test should run on. */
#define MIN_NTASKS 4

/* The name of this test. */
#define TEST_NAME "test_darray_multivar"

/* Number of processors that will do IO. */
#define NUM_IO_PROCS 1

/* Number of computational components to create. */
#define COMPONENT_COUNT 1

/* The number of dimensions in the example data. In this test, we
 * are using three-dimensional data. */
#define NDIM 3

/* But sometimes we need arrays of the non-record dimensions. */
#define NDIM2 2

/* The length of our sample data along each dimension. */
#define X_DIM_LEN 4
#define Y_DIM_LEN 4

/* The number of timesteps of data to write. */
#define NUM_TIMESTEPS 2

/* The name of the variables in the netCDF output files. */
#define VAR_NAME_1 "STICKS"
#define VAR_NAME_2 "NIX"
#define VAR_NAME_3 "HICK"
#define VAR_NAME_4 "PIX"

/* Number of variables in the test file. */
#define NUM_VAR 4

/* The dimension names. */
char dim_name[NDIM][PIO_MAX_NAME + 1] = {"timestep", "x", "y"};

/* The var names. */
char var_name[NUM_VAR][PIO_MAX_NAME + 1] = {"STICKS", "NIX", "HICK", "PIX"};

/* Length of the dimensions in the sample data. */
int dim_len[NDIM] = {NC_UNLIMITED, X_DIM_LEN, Y_DIM_LEN};

/** 
 * Test the darray functionality. Create a netCDF file with 3
 * dimensions and 4 variables, and use darray to write to one of them.
 *
 * @param iosysid the IO system ID.
 * @param ioid the ID of the decomposition.
 * @param num_flavors the number of IOTYPES available in this build.
 * @param flavor array of available iotypes.
 * @param my_rank rank of this task.
 * @param pio_type the type of the data.
 * @param test_comm the communicator that is running this test.
 * @param rearranger the rearranger in use for this test.
 * @param use_fill 1 if fill mode should be set.
 * @param use_default 1 if default fill values should be used
 * (ignored if use_fill is 0).
 * @returns 0 for success, error code otherwise.
*/
int test_3_empty(int iosysid, int ioid, int num_flavors, int *flavor, int my_rank,
                 int pio_type, MPI_Comm test_comm, int rearranger, int use_fill,
                 int use_default)
{
    char filename[PIO_MAX_NAME + 1]; /* Name for the output files. */
    int dimids[NDIM];     /* The dimension IDs. */
    int ncid;             /* The ncid of the netCDF file. */
    int ncid2;            /* The ncid of the re-opened netCDF file. */
    int varid[NUM_VAR];   /* The IDs of the netCDF varables. */
    PIO_Offset arraylen = 4;
    void *fillvalue;
    void *test_data;
    void *test_data_in;
    int fillvalue_int = NC_FILL_INT;
    int custom_fillvalue_int = -TEST_VAL_42;
    int test_data_int[arraylen];
    int test_data_int_in[arraylen];
    float fillvalue_float = NC_FILL_FLOAT;
    float custom_fillvalue_float = -TEST_VAL_42;
    float test_data_float[arraylen];
    float test_data_float_in[arraylen];
    double fillvalue_double = NC_FILL_DOUBLE;
    double custom_fillvalue_double = (-TEST_VAL_42 * 100);
    double test_data_double[arraylen];
    double test_data_double_in[arraylen];
    int ret;       /* Return code. */

    /* Initialize some data. */
    for (int f = 0; f < arraylen; f++)
    {
        test_data_int[f] = my_rank * 10 + f;
        test_data_float[f] = my_rank * 10 + f + 0.5;
        test_data_double[f] = my_rank * 100000 + f + 0.5;
    }

    /* Select the fill value and data. */
    switch (pio_type)
    {
    case PIO_INT:
        fillvalue = use_default ? &fillvalue_int : &custom_fillvalue_int;
        test_data = test_data_int;
        test_data_in = test_data_int_in;
        break;
    case PIO_FLOAT:
        fillvalue = use_default ? &fillvalue_float : &custom_fillvalue_float;
        test_data = test_data_float;
        test_data_in = test_data_float_in;
        break;
    case PIO_DOUBLE:
        fillvalue = use_default ? &fillvalue_double : &custom_fillvalue_double;
        test_data = test_data_double;
        test_data_in = test_data_double_in;
        break;
    default:
        ERR(ERR_WRONG);
    }

    /* Try in pnetcdf only. This code demonstrates that pnetcdf fill
     * values do work. But not for PIO, so we have a bug somewhere. */
    {
#ifdef _PNETCDF
        int ncid;
        int varid;
        int dimid;
        char test_filename[] = "pnetcdf_test.nc";
        int ret;

        if ((ret = ncmpi_create(test_comm, test_filename, NC_CLOBBER, MPI_INFO_NULL, &ncid)))
            return ret;
        if ((ret = ncmpi_set_fill(ncid, NC_FILL, NULL)))
            return ret;            
        if ((ret = ncmpi_def_dim(ncid, "dim_name", 5, &dimid)))
            return ret;
        if ((ret = ncmpi_def_var(ncid, "dim_name", NC_INT, 1, &dimid, &varid)))
            return ret;
        if ((ret = ncmpi_enddef(ncid)))
            return ret;
        if ((ret = ncmpi_close(ncid)))
            return ret;

        /* Reopen and check. */
        if ((ret = ncmpi_open(test_comm, test_filename, NC_NOWRITE, MPI_INFO_NULL, &ncid)))
            return ret;
        int datum;
        MPI_Offset start[1] = {0};
        ret = ncmpi_get_var1_int(ncid, varid, start, &datum);
        printf("datum ret = %d\n", ret);

        /* Not sure why this doesn't work. */
        /* if ((ret = ncmpi_get_var1_int(ncid, varid, start, &datum))) */
        /*     return ret; */
        printf("datum = %d\n", datum);
        if ((ret = ncmpi_close(ncid)))
            return ret;
#endif /* _PNETCDF */
    }

    /* Use PIO to create the example file in each of the four
     * available ways. */
    for (int fmt = 0; fmt < num_flavors; fmt++) 
    {
        /* Create the filename. */
        sprintf(filename, "data_%s_iotype_%d_pio_type_%d_use_fill_%d_default_fill_%d.nc",
                TEST_NAME, flavor[fmt], pio_type, use_fill, use_default);

        /* Create the netCDF output file. */
        printf("rank: %d Creating sample file %s with format %d type %d\n", my_rank, filename,
               flavor[fmt], pio_type);
        if ((ret = PIOc_createfile(iosysid, &ncid, &flavor[fmt], filename, PIO_CLOBBER)))
            ERR(ret);

        /* Turn on fill mode if desired. */
        if (use_fill)
            if ((ret = PIOc_set_fill(ncid, NC_FILL, NULL)))
                ERR(ret);

        /* Define netCDF dimensions and variable. */
        printf("%d Defining netCDF metadata...\n", my_rank);
        for (int d = 0; d < NDIM; d++)
            if ((ret = PIOc_def_dim(ncid, dim_name[d], (PIO_Offset)dim_len[d], &dimids[d])))
                ERR(ret);

        /* Define the variables. */
        for (int v = 0; v < NUM_VAR; v++)
        {
            if ((ret = PIOc_def_var(ncid, var_name[v], pio_type, NDIM, dimids, &varid[v])))
                ERR(ret);
            if (use_fill && !use_default)
                if ((ret = PIOc_def_var_fill(ncid, varid[v], NC_FILL, fillvalue)))
                    ERR(ret);
        }

        /* End define mode. */
        if ((ret = PIOc_enddef(ncid)))
            ERR(ret);

        /* Set the value of the record dimension. */
        if ((ret = PIOc_setframe(ncid, varid[0], 0)))
            ERR(ret);

        /* Write the data. */
        if ((ret = PIOc_write_darray(ncid, varid[0], ioid, arraylen, test_data, fillvalue)))
            ERR(ret);

        /* Close the netCDF file. */
        if ((ret = PIOc_closefile(ncid)))
            ERR(ret);

        /* Reopen the file. */
        if ((ret = PIOc_openfile(iosysid, &ncid2, &flavor[fmt], filename, PIO_NOWRITE)))
            ERR(ret);

        /* Read the data. */
        if ((ret = PIOc_read_darray(ncid2, varid[0], ioid, arraylen, test_data_in)))
            ERR(ret);

        /* Check the results. */
        for (int f = 0; f < arraylen; f++)
        {
            switch (pio_type)
            {
            case PIO_INT:
                if (test_data_int_in[f] != test_data_int[f])
                    return ERR_WRONG;
                break;
            case PIO_FLOAT:
                if (test_data_float_in[f] != test_data_float[f])
                    return ERR_WRONG;
                break;
            case PIO_DOUBLE:
                if (test_data_double_in[f] != test_data_double[f])
                    return ERR_WRONG;
                break;
            default:
                ERR(ERR_WRONG);
            }
        }

        /* If fill mode is in use the other vars should have fill values. */
        if (use_fill && flavor[fmt] != PIO_IOTYPE_PNETCDF)
        {
            /* Read the data. */
            if ((ret = PIOc_read_darray(ncid2, varid[1], ioid, arraylen, test_data_in)))
                ERR(ret);
            
            /* Check the results. */
            for (int f = 0; f < arraylen; f++)
            {
                switch (pio_type)
                {
                case PIO_INT:
                    if (test_data_int_in[f] != (use_default ? NC_FILL_INT : custom_fillvalue_int))
                        return ERR_WRONG;
                    break;
                case PIO_FLOAT:
                    if (test_data_float_in[f] != (use_default ? NC_FILL_FLOAT : custom_fillvalue_float))
                        return ERR_WRONG;
                    break;
                case PIO_DOUBLE:
                    if (test_data_double_in[f] != (use_default ? NC_FILL_DOUBLE : custom_fillvalue_double))
                        return ERR_WRONG;
                    break;
                default:
                    ERR(ERR_WRONG);
                }
            }
        }

        /* Close the netCDF file. */
        printf("%d Closing the sample data file...\n", my_rank);
        if ((ret = PIOc_closefile(ncid2)))
            ERR(ret);
    }
    return PIO_NOERR;
}

/**
 * Run all the tests. 
 *
 * @param iosysid the IO system ID.
 * @param num_flavors number of available iotypes in the build.
 * @param flavor pointer to array of the available iotypes.
 * @param my_rank rank of this task.
 * @param test_comm the communicator the test is running on.
 * @param rearranger the rearranger in use in this test.
 * @returns 0 for success, error code otherwise. 
 */
int test_all_darray(int iosysid, int num_flavors, int *flavor, int my_rank,
                    MPI_Comm test_comm, int rearranger)
{
#define NUM_FILL_TESTS 3
#define NUM_TYPES_TO_TEST 3
    int pio_type[NUM_TYPES_TO_TEST] = {PIO_INT, PIO_FLOAT, PIO_DOUBLE};
    int ioid;
    int dim_len_2d[NDIM2] = {X_DIM_LEN, Y_DIM_LEN};
    int ret; /* Return code. */

    for (int t = 0; t < NUM_TYPES_TO_TEST; t++)
    {        
        int use_fill = 0;
        int use_default = 0;
        
        /* Decompose the data over the tasks. */
        if ((ret = create_decomposition_2d(TARGET_NTASKS, my_rank, iosysid, dim_len_2d,
                                           &ioid, pio_type[t])))
        return ret;

        /* Run the different combinations of use_fill and use_default. */
        for (int f = 0; f < NUM_FILL_TESTS; f++)
        {
            /* Set flags for this test case. */
            if (f == 1)
                use_fill++;
            if (f == 2)
                use_default++;
            
            /* Run a simple darray test. */
            if ((ret = test_3_empty(iosysid, ioid, num_flavors, flavor, my_rank, pio_type[t],
                                    test_comm, rearranger, use_fill, use_default)))
                return ret;
        }
        
        /* Free the PIO decomposition. */
        if ((ret = PIOc_freedecomp(iosysid, ioid)))
            ERR(ret);
    }

    return PIO_NOERR;
}

/* Run tests for darray functions. */
int main(int argc, char **argv)
{
#define NUM_REARRANGERS 2
    int rearranger[NUM_REARRANGERS] = {PIO_REARR_BOX, PIO_REARR_SUBSET};
    int my_rank;
    int ntasks;
    int num_flavors;         /* Number of PIO netCDF flavors in this build. */
    int flavor[NUM_FLAVORS]; /* iotypes for the supported netCDF IO flavors. */
    MPI_Comm test_comm;      /* A communicator for this test. */
    int ret;                 /* Return code. */

    /* Initialize test. */
    if ((ret = pio_test_init2(argc, argv, &my_rank, &ntasks, MIN_NTASKS, MIN_NTASKS,
                              3, &test_comm)))
        ERR(ERR_INIT);

    if ((ret = PIOc_set_iosystem_error_handling(PIO_DEFAULT, PIO_RETURN_ERROR, NULL)))
        return ret;

    /* Only do something on max_ntasks tasks. */
    if (my_rank < TARGET_NTASKS)
    {
        int iosysid;              /* The ID for the parallel I/O system. */
        int ioproc_stride = 1;    /* Stride in the mpi rank between io tasks. */
        int ioproc_start = 0;     /* Zero based rank of first processor to be used for I/O. */
        int ret;                  /* Return code. */
        
        /* Figure out iotypes. */
        if ((ret = get_iotypes(&num_flavors, flavor)))
            ERR(ret);
        printf("Runnings tests for %d flavors\n", num_flavors);

        /* Test for both arrangers. */
        for (int r = 0; r < NUM_REARRANGERS; r++)
        {

            /* Initialize the PIO IO system. This specifies how
             * many and which processors are involved in I/O. */
            if ((ret = PIOc_Init_Intracomm(test_comm, TARGET_NTASKS, ioproc_stride,
                                           ioproc_start, rearranger[r], &iosysid)))
                return ret;
            
            /* Run tests. */
            printf("%d Running tests...\n", my_rank);
            if ((ret = test_all_darray(iosysid, num_flavors, flavor, my_rank, test_comm,
                                       rearranger[r])))
                return ret;
            
            /* Finalize PIO system. */
            if ((ret = PIOc_finalize(iosysid)))
                return ret;
        }

    } /* endif my_rank < TARGET_NTASKS */

    /* Finalize the MPI library. */
    printf("%d %s Finalizing...\n", my_rank, TEST_NAME);
    if ((ret = pio_test_finalize(&test_comm)))
        return ret;

    printf("%d %s SUCCESS!!\n", my_rank, TEST_NAME);
    return 0;
}
