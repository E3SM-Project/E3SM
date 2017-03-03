/*
 * Tests for PIO distributed arrays. This test uses 1 dimension,
 * everything very simple. ;-)
 *
 * Ed Hartnett, 2/27/17
 */
#include <pio.h>
#include <pio_internal.h>
#include <pio_tests.h>

/* The number of tasks this test should run on. */
#define TARGET_NTASKS 4

/* The minimum number of tasks this test should run on. */
#define MIN_NTASKS 4

/* The name of this test. */
#define TEST_NAME "test_darray_1d"

/* Number of processors that will do IO. */
#define NUM_IO_PROCS 1

/* Number of computational components to create. */
#define COMPONENT_COUNT 1

/* The number of dimensions in the example data. In this test, we
 * are using three-dimensional data. */
#define NDIM 1

/* The length of our sample data along the dimension. */
#define DIM_LEN 8

/* This is the length of the map for each task. */
#define EXPECTED_MAPLEN 2

/* The number of timesteps of data to write. */
#define NUM_TIMESTEPS 2

/* The name of the variable in the netCDF output files. */
#define VAR_NAME "RedShirtSurvival"

/* The dimension names. */
#define DIM_NAME "episode"
#define DIM_NAME_2 "phaser_draws"

/* Create a 1D decomposition.
 *
 * @param ntasks the number of available tasks
 * @param my_rank rank of this task.
 * @param iosysid the IO system ID.
 * @param dim_len an array of length 3 with the dimension sizes.
 * @param ioid a pointer that gets the ID of this decomposition.
 * @param pio_type the type that will be used for basetype.
 * @returns 0 for success, error code otherwise.
 **/
int create_decomposition_1d(int ntasks, int my_rank, int iosysid, int pio_type, int *ioid)
{
    PIO_Offset elements_per_pe;     /* Array elements per processing unit. */
    int dim_len_1d[NDIM] = {DIM_LEN};
    int ret;

    /* How many data elements per task? In this example we will end up
     * with 2. */
    elements_per_pe = DIM_LEN / ntasks;

    PIO_Offset compdof[elements_per_pe];

    /* Don't forget to add 1! */
    compdof[0] = my_rank + 1;

    /* This means fill value will be used here. */
    compdof[1] = 0;

    /* Create the PIO decomposition for this test. */
    if ((ret = PIOc_InitDecomp(iosysid, pio_type, NDIM, dim_len_1d, elements_per_pe,
                               compdof, ioid, NULL, NULL, NULL)))
        ERR(ret);

    printf("%d decomposition initialized.\n", my_rank);

    return 0;
}

/**
 * Test fill values and darrays.
 *
 * @param iosysid the IO system ID.
 * @param ioid the ID of the decomposition.
 * @param pio_type the type of the data.
 * @param num_flavors the number of IOTYPES available in this build.
 * @param flavor array of available iotypes.
 * @param my_rank rank of this task.
 * @param test_comm the MPI communicator running the test.
 * @returns 0 for success, error code otherwise.
*/
int test_darray_fill(int iosysid, int ioid, int pio_type, int num_flavors, int *flavor,
                     int my_rank, MPI_Comm test_comm)
{
#define NUM_FILLVALUE_PRESENT_TESTS 2
    char filename[PIO_MAX_NAME + 1]; /* Name for the output files. */
    int dimid;     /* The dimension ID. */
    int ncid;      /* The ncid of the netCDF file. */
    int varid;     /* The ID of the netCDF varable. */
    PIO_Offset arraylen = 2;
    void *test_data;
    void *fillvalue;
    void *test_data_in;
    void *expected_in;
    PIO_Offset type_size;             /* Size of the data type. */
    float my_float_rank = my_rank;    /* my_rank in a float. */
    double my_double_rank = my_rank;  /* my_rank in a double. */
    int int_fill = NC_FILL_INT;
    float float_fill = NC_FILL_FLOAT;
    double double_fill = NC_FILL_DOUBLE;
    void *bufr;
    int ret;                        /* Return code. */

    /* Use PIO to create the example file in each of the four
     * available ways. */
    for (int fmt = 0; fmt < num_flavors; fmt++)
    {
        for (int with_fillvalue = 0; with_fillvalue < NUM_FILLVALUE_PRESENT_TESTS; with_fillvalue++)
        {
            /* Create the filename. */
            sprintf(filename, "data_%s_iotype_%d_pio_type_%d_with_fillvalue_%d.nc", TEST_NAME, flavor[fmt],
                    pio_type, with_fillvalue);

            /* Create the netCDF output file. */
            printf("rank: %d Creating sample file %s with format %d...\n", my_rank, filename,
                   flavor[fmt]);
            if ((ret = PIOc_createfile(iosysid, &ncid, &flavor[fmt], filename, PIO_CLOBBER)))
                ERR(ret);

            /* Turn on fill mode. */
            if ((ret = PIOc_set_fill(ncid, NC_FILL, NULL)))
                ERR(ret);

            /* Define netCDF dimensions and variable. */
            if ((ret = PIOc_def_dim(ncid, DIM_NAME, DIM_LEN, &dimid)))
                ERR(ret);

            /* Define a variable. */
            if ((ret = PIOc_def_var(ncid, VAR_NAME, pio_type, NDIM, &dimid, &varid)))
                ERR(ret);

            /* End define mode. */
            if ((ret = PIOc_enddef(ncid)))
                ERR(ret);

            /* Get the size of the type. */
            if ((ret = PIOc_inq_type(ncid, pio_type, NULL, &type_size)))
                return ret;

            /* Initialize some data. */
            int int_test_data[2] = {my_rank, my_rank};
            float float_test_data[2] = {my_rank, my_rank};
            double double_test_data[2] = {my_rank, my_rank};
            switch (pio_type)
            {
            case PIO_INT:
                test_data = int_test_data;
                fillvalue = with_fillvalue ? &int_fill : NULL;
                expected_in = &my_rank;
                break;
            case PIO_FLOAT:
                test_data = float_test_data;
                fillvalue = with_fillvalue ? &float_fill : NULL;
                expected_in = &my_float_rank;
                break;
            case PIO_DOUBLE:
                test_data = double_test_data;
                fillvalue = with_fillvalue ? &double_fill : NULL;
                expected_in = &my_double_rank;
                break;
            default:
                return ERR_WRONG;
            }

            /* Write the data. Our test_data contains only one real value
             * (instead of 2, as indicated by arraylen), but due to the
             * decomposition, only the first value is used in the
             * output. */
            if ((ret = PIOc_write_darray(ncid, varid, ioid, arraylen, test_data,
                                         fillvalue)))
                ERR(ret);

            /* Close the netCDF file. */
            if ((ret = PIOc_closefile(ncid)))
                ERR(ret);

            /* Reopen the file. */
            if ((ret = PIOc_openfile(iosysid, &ncid, &flavor[fmt], filename, PIO_NOWRITE)))
                ERR(ret);

            /* Allocate space for data. */
            if (!(test_data_in = malloc(type_size * arraylen)))
                ERR(PIO_ENOMEM);

            /* Read the data. */
            if ((ret = PIOc_read_darray(ncid, varid, ioid, arraylen, test_data_in)))
                ERR(ret);

            /* Check the (first) result. */
            if (memcmp(test_data_in, expected_in, type_size))
                return ERR_WRONG;

            /* Free resources. */
            free(test_data_in);

            /* Get a buffer big enough to hold the global array. */
            if (!(bufr = malloc(DIM_LEN * type_size)))
                return PIO_ENOMEM;

            /* Get the whole array with good old get_var(). */
            if ((ret = PIOc_get_var(ncid, varid, bufr)))
                return ret;

            /* Check the results. The first four values are 0, 1, 2, 3,
             * and the rest are the default fill value of the type. */
            for (int e = 0; e < DIM_LEN; e++)
            {
                switch (pio_type)
                {
                case PIO_INT:
                    if (((int *)bufr)[e] != (e < 4 ? e : NC_FILL_INT))
                        return ERR_WRONG;
                    break;
                case PIO_FLOAT:
                    if (((float *)bufr)[e] != (e < 4 ? e : NC_FILL_FLOAT))
                        return ERR_WRONG;
                    break;
                case PIO_DOUBLE:
                    if (((double *)bufr)[e] != (e < 4 ? e : NC_FILL_DOUBLE))
                        return ERR_WRONG;
                    break;
                default:
                    return ERR_WRONG;
                }
            }

            /* Release buffer. */
            free(bufr);

            /* Close the netCDF file. */
            printf("%d Closing the sample data file...\n", my_rank);
            if ((ret = PIOc_closefile(ncid)))
                ERR(ret);
        } /* with_fillvalue */
    } /* next iotype */

    return PIO_NOERR;
}

/**
 * Test fill values and darrays with an unlimited dim.
 *
 * @param iosysid the IO system ID.
 * @param ioid the ID of the decomposition.
 * @param pio_type the type of the data.
 * @param num_flavors the number of IOTYPES available in this build.
 * @param flavor array of available iotypes.
 * @param my_rank rank of this task.
 * @param test_comm the MPI communicator running the test.
 * @returns 0 for success, error code otherwise.
*/
int test_darray_fill_unlim(int iosysid, int ioid, int pio_type, int num_flavors,
                           int *flavor, int my_rank, MPI_Comm test_comm)
{
#define NDIM2 2
    char filename[PIO_MAX_NAME + 1]; /* Name for the output files. */
    int dimid[NDIM2];     /* The dimension ID. */
    int ncid;      /* The ncid of the netCDF file. */
    int varid;     /* The ID of the netCDF varable. */
    PIO_Offset arraylen = 2;
    void *test_data;
    void *fillvalue;
    void *test_data_in;
    void *expected_in;
    PIO_Offset type_size;             /* Size of the data type. */
    float my_float_rank = my_rank;    /* my_rank in a float. */
    double my_double_rank = my_rank;  /* my_rank in a double. */
    int int_fill = NC_FILL_INT;
    float float_fill = NC_FILL_FLOAT;
    double double_fill = NC_FILL_DOUBLE;
    void *bufr;
    int ret;                        /* Return code. */

    /* Use PIO to create the example file in each of the four
     * available ways. */
    for (int fmt = 0; fmt < num_flavors; fmt++)
    {
        /* Create the filename. */
        sprintf(filename, "data_%s_iotype_%d_pio_type_%d_unlim.nc", TEST_NAME, flavor[fmt],
                pio_type);

        /* Create the netCDF output file. */
        printf("rank: %d Creating sample file %s with format %d...\n", my_rank, filename,
               flavor[fmt]);
        if ((ret = PIOc_createfile(iosysid, &ncid, &flavor[fmt], filename, PIO_CLOBBER)))
            ERR(ret);

        /* Turn on fill mode. */
        if ((ret = PIOc_set_fill(ncid, NC_FILL, NULL)))
            ERR(ret);

        /* Define netCDF dimensions. */
        if ((ret = PIOc_def_dim(ncid, DIM_NAME, NC_UNLIMITED, &dimid[0])))
            ERR(ret);
        if ((ret = PIOc_def_dim(ncid, DIM_NAME_2, DIM_LEN, &dimid[1])))
            ERR(ret);

        /* Define a variable. */
        if ((ret = PIOc_def_var(ncid, VAR_NAME, pio_type, NDIM2, dimid, &varid)))
            ERR(ret);

        /* End define mode. */
        if ((ret = PIOc_enddef(ncid)))
            ERR(ret);

        /* Get the size of the type. */
        if ((ret = PIOc_inq_type(ncid, pio_type, NULL, &type_size)))
            return ret;

        /* Initialize some data. */
        int int_test_data[2] = {my_rank, my_rank};
        float float_test_data[2] = {my_rank, my_rank};
        double double_test_data[2] = {my_rank, my_rank};
        switch (pio_type)
        {
        case PIO_INT:
            test_data = int_test_data;
            fillvalue = &int_fill;
            expected_in = &my_rank;
            break;
        case PIO_FLOAT:
            test_data = float_test_data;
            fillvalue = &float_fill;
            expected_in = &my_float_rank;
            break;
        case PIO_DOUBLE:
            test_data = double_test_data;
            fillvalue = &double_fill;
            expected_in = &my_double_rank;
            break;
        default:
            return ERR_WRONG;
        }

        /* Set the record number for the unlimited dimension. */
        if ((ret = PIOc_setframe(ncid, varid, 0)))
            ERR(ret);

        /* Write the data. Our test_data contains only one real value
         * (instead of 2, as indicated by arraylen), but due to the
         * decomposition, only the first value is used in the
         * output. */
        if ((ret = PIOc_write_darray(ncid, varid, ioid, arraylen, test_data, fillvalue)))
            ERR(ret);

        /* Set the record number for the unlimited dimension. */
        if ((ret = PIOc_setframe(ncid, varid, 1)))
            ERR(ret);

        /* Write the data. Our test_data contains only one real value
         * (instead of 2, as indicated by arraylen), but due to the
         * decomposition, only the first value is used in the
         * output. */
        if ((ret = PIOc_write_darray(ncid, varid, ioid, arraylen, test_data, fillvalue)))
            ERR(ret);

        /* Close the netCDF file. */
        if ((ret = PIOc_closefile(ncid)))
            ERR(ret);

        /* Reopen the file. */
        if ((ret = PIOc_openfile(iosysid, &ncid, &flavor[fmt], filename, PIO_NOWRITE)))
            ERR(ret);

        /* Allocate space for data. */
        if (!(test_data_in = malloc(type_size * arraylen)))
            ERR(PIO_ENOMEM);

        /* Read the data. */
        if ((ret = PIOc_read_darray(ncid, varid, ioid, arraylen, test_data_in)))
            ERR(ret);

        /* Check the (first) result. */
        if (memcmp(test_data_in, expected_in, type_size))
            return ERR_WRONG;

        /* Free resources. */
        free(test_data_in);

        /* Get a buffer big enough to hold the global array. */
        if (!(bufr = malloc(DIM_LEN * type_size * 2)))
            return PIO_ENOMEM;

        /* Get the whole array with good old get_var(). */
        if ((ret = PIOc_get_var(ncid, varid, bufr)))
            return ret;

        /* Check the results. The first four values in each record are
         * 0, 1, 2, 3, and the rest are the default fill value of the
         * type. There are two records. */
        for (int e = 0; e < DIM_LEN * 2; e++)
        {
            switch (pio_type)
            {
            case PIO_INT:
                if (((int *)bufr)[e] != (e % 8 < 4 ? e % 8 : NC_FILL_INT))
                    return ERR_WRONG;
                break;
            case PIO_FLOAT:
                if (((float *)bufr)[e] != (e % 8 < 4 ? e % 8 : NC_FILL_FLOAT))
                    return ERR_WRONG;
                break;
            case PIO_DOUBLE:
                if (((double *)bufr)[e] != (e % 8 < 4 ? e % 8 : NC_FILL_DOUBLE))
                    return ERR_WRONG;
                break;
            default:
                return ERR_WRONG;
            }
        }

        /* Release buffer. */
        free(bufr);

        /* Close the netCDF file. */
        printf("%d Closing the sample data file...\n", my_rank);
        if ((ret = PIOc_closefile(ncid)))
            ERR(ret);
    } /* next iotype */

    return PIO_NOERR;
}

/**
 * Test the decomp read/write functionality.
 *
 * @param iosysid the IO system ID.
 * @param ioid the ID of the decomposition.
 * @param num_flavors the number of IOTYPES available in this build.
 * @param flavor array of available iotypes.
 * @param my_rank rank of this task.
 * @param pio_type the type involved in this decompositon.
 * @param rearranger the rearranger in use.
 * @param test_comm the MPI communicator for this test.
 * @returns 0 for success, error code otherwise.
*/
int test_decomp_read_write(int iosysid, int ioid, int num_flavors, int *flavor, int my_rank,
                           int pio_type, int rearranger, MPI_Comm test_comm)
{
    char filename[PIO_MAX_NAME + 1]; /* Name for the output files. */
    int ioid2;             /* ID for decomposition we will create from file. */
    char title_in[PIO_MAX_NAME + 1];   /* Optional title. */
    char history_in[PIO_MAX_NAME + 1]; /* Optional history. */
    int fortran_order_in; /* Indicates fortran vs. c order. */
    int ret;              /* Return code. */

    /* Use PIO to create the decomp file in each of the four
     * available ways. */
    for (int fmt = 0; fmt < num_flavors; fmt++)
    {
        /* Create the filename. */
        sprintf(filename, "decomp_%s_iotype_%d.nc", TEST_NAME, flavor[fmt]);

        printf("writing decomp file %s\n", filename);
        if ((ret = PIOc_write_nc_decomp(iosysid, filename, 0, ioid, test_comm, NULL,
                                        NULL, 0)))
            return ret;

        /* Read the data. */
        printf("reading decomp file %s\n", filename);
        if ((ret = PIOc_read_nc_decomp(iosysid, filename, &ioid2, test_comm, pio_type,
                                       title_in, history_in, &fortran_order_in)))
            return ret;

        /* Check the results. */
        {
            iosystem_desc_t *ios;
            io_desc_t *iodesc;
            int expected_basetype;

            switch (pio_type)
            {
            case PIO_INT:
                expected_basetype = MPI_INT;
                break;
            case PIO_FLOAT:
                expected_basetype = MPI_FLOAT;
                break;
            case PIO_DOUBLE:
                expected_basetype = MPI_DOUBLE;
                break;
            default:
                return ERR_WRONG;
            }

            /* Get the IO system info. */
            if (!(ios = pio_get_iosystem_from_id(iosysid)))
                return pio_err(NULL, NULL, PIO_EBADID, __FILE__, __LINE__);

            /* Get the IO desc, which describes the decomposition. */
            if (!(iodesc = pio_get_iodesc_from_id(ioid2)))
                return pio_err(ios, NULL, PIO_EBADID, __FILE__, __LINE__);
            if (iodesc->ioid != ioid2 || iodesc->maplen != EXPECTED_MAPLEN || iodesc->ndims != NDIM)
                return ERR_WRONG;
            /* if (iodesc->nrecvs != 1) */
            /*     return ERR_WRONG; */
            /* if (iodesc->num_aiotasks != TARGET_NTASKS) */
            /*     return ERR_WRONG; */
            printf("iodesc->nrecvs = %d iodesc->num_aiotasks = %d\n", iodesc->nrecvs, iodesc->num_aiotasks);
            if (iodesc->ndof != EXPECTED_MAPLEN)
                return ERR_WRONG;
            if (iodesc->rearranger != rearranger || iodesc->maxregions != 1)
                return ERR_WRONG;
            if (!iodesc->needsfill || iodesc->basetype != expected_basetype)
                return ERR_WRONG;
            /* Don't forget to add 1!! */
            if (iodesc->map[0] != my_rank + 1 || iodesc->map[1] != 0)
                return ERR_WRONG;
            if (iodesc->dimlen[0] != DIM_LEN)
                return ERR_WRONG;
        }

        /* Free the PIO decomposition. */
        if ((ret = PIOc_freedecomp(iosysid, ioid2)))
            ERR(ret);
    }
    return PIO_NOERR;
}

/* Run tests for darray functions. */
int main(int argc, char **argv)
{
#define NUM_REARRANGERS_TO_TEST 2
    int rearranger[NUM_REARRANGERS_TO_TEST] = {PIO_REARR_BOX, PIO_REARR_SUBSET};
#define NUM_TYPES_TO_TEST 3
    int test_type[NUM_TYPES_TO_TEST] = {PIO_INT, PIO_FLOAT, PIO_DOUBLE};
    int my_rank;
    int ntasks;
    int num_flavors; /* Number of PIO netCDF flavors in this build. */
    int flavor[NUM_FLAVORS]; /* iotypes for the supported netCDF IO flavors. */
    MPI_Comm test_comm; /* A communicator for this test. */
    int ret;         /* Return code. */

    /* Initialize test. */
    if ((ret = pio_test_init2(argc, argv, &my_rank, &ntasks, MIN_NTASKS,
                              MIN_NTASKS, 3, &test_comm)))
        ERR(ERR_INIT);

    if ((ret = PIOc_set_iosystem_error_handling(PIO_DEFAULT, PIO_RETURN_ERROR, NULL)))
        return ret;

    /* Only do something on max_ntasks tasks. */
    if (my_rank < TARGET_NTASKS)
    {
        int iosysid;  /* The ID for the parallel I/O system. */
        int ioid;     /* Decomposition ID. */
        int ioproc_stride = 1;    /* Stride in the mpi rank between io tasks. */
        int ioproc_start = 0;     /* Rank of first processor to be used for I/O. */
        int ret;      /* Return code. */

        /* Figure out iotypes. */
        if ((ret = get_iotypes(&num_flavors, flavor)))
            ERR(ret);
        printf("Runnings tests for %d flavors\n", num_flavors);

        for (int r = 0; r < NUM_REARRANGERS_TO_TEST; r++)
        {
            /* Initialize the PIO IO system. This specifies how many and
             * which processors are involved in I/O. */
            if ((ret = PIOc_Init_Intracomm(test_comm, TARGET_NTASKS, ioproc_stride,
                                           ioproc_start, rearranger[r], &iosysid)))
                return ret;

            /* Run tests for each data type. */
            for (int t = 0; t < NUM_TYPES_TO_TEST; t++)
            {
                /* Decompose the data over the tasks. */
                if ((ret = create_decomposition_1d(TARGET_NTASKS, my_rank, iosysid, test_type[t],
                                                   &ioid)))
                    return ret;

                /* Test decomposition read/write. */
                if ((ret = test_decomp_read_write(iosysid, ioid, num_flavors, flavor, my_rank,
                                                  test_type[t], rearranger[r], test_comm)))
                    return ret;

                /* Run tests. */
                if ((ret = test_darray_fill(iosysid, ioid, test_type[t], num_flavors, flavor,
                                            my_rank, test_comm)))
                    return ret;

                /* Run tests. */
                if ((ret = test_darray_fill_unlim(iosysid, ioid, test_type[t], num_flavors,
                                                  flavor, my_rank, test_comm)))
                    return ret;

                /* Free the PIO decomposition. */
                if ((ret = PIOc_freedecomp(iosysid, ioid)))
                    ERR(ret);
            }

            /* Finalize PIO system. */
            if ((ret = PIOc_finalize(iosysid)))
                return ret;
        } /* next rearranger */

    } /* endif my_rank < TARGET_NTASKS */

    /* Finalize the MPI library. */
    printf("%d %s Finalizing...\n", my_rank, TEST_NAME);
    if ((ret = pio_test_finalize(&test_comm)))
        return ret;

    printf("%d %s SUCCESS!!\n", my_rank, TEST_NAME);
    return 0;
}
