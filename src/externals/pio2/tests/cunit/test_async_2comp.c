/**
 * @file Tests for PIOc_Intercomm. This tests basic asynch I/O capability.
 * @author Ed Hartnett
 *
 * To run with valgrind, use this command:
 * <pre>mpiexec -n 4 valgrind -v --leak-check=full --suppressions=../../../tests/unit/valsupp_test.supp
 * --error-exitcode=99 --track-origins=yes ./test_intercomm3</pre>
 *
 */
#include <pio.h>
#include <pio_tests.h>

/* Number of processors that will do IO. */
#define NUM_IO_PROCS 2

/* Number of computational components to create. */
#define COMPONENT_COUNT 2

/* The number of tasks this test should run on. */
#define TARGET_NTASKS 4

/* The name of this test. */
#define TEST_NAME "test_intercomm3"

/** Run Tests for Init_Intercomm
 *
 */
int
main(int argc, char **argv)
{
    int verbose = 1;

    /* Zero-based rank of processor. */
    int my_rank;

    /* Number of processors involved in current execution. */
    int ntasks;

    /* Different output flavors. */
    int flavor[NUM_FLAVORS];

    int num_flavors;

    /* The ID for the parallel I/O system. */
    int iosysid[COMPONENT_COUNT];

    /* The ncid of the netCDF file. */
    int ncid;

    /* The ID of the netCDF varable. */
    int varid;

    /* Return code. */
    int ret;

    /* Index for loops. */
    int fmt, d, d1, i;

    /* Initialize test. */
    if ((ret = pio_test_init(argc, argv, &my_rank, &ntasks, TARGET_NTASKS)))
        ERR(ERR_INIT);

    /* Figure out iotypes. */
    if ((ret = get_iotypes(&num_flavors, flavor)))
        ERR(ret);

    /* How many processors will be used for our IO and 2 computation components. */
    int num_procs[COMPONENT_COUNT + 1] = {2, 1, 1};

    /* Is the current process a computation task? */
    int comp_task = my_rank < 2 ? 0 : 1;

    /* Index of computation task in iosysid array. Varies by rank and
     * does not apply to IO component processes. */
    int my_comp_idx = comp_task ? my_rank - 2 : -1;

    /* Initialize the IO system. */
    if ((ret = PIOc_Init_Async(MPI_COMM_WORLD, NUM_IO_PROCS, NULL, COMPONENT_COUNT,
                               num_procs, NULL, iosysid)))
        ERR(ERR_AWFUL);

    /* All the netCDF calls are only executed on the computation
     * tasks. The IO tasks have not returned from PIOc_Init_Intercomm,
     * and when the do, they should go straight to finalize. */
    if (comp_task)
    {
        for (int flv = 0; flv < num_flavors; flv++)
        {
            char filename[NC_MAX_NAME + 1];

            /* Create a filename. */
            int sample = 1;
            sprintf(filename, "%s_%s_%d_%d.nc", TEST_NAME, flavor_name(flv), sample, my_comp_idx);

            /* Create sample file 1. */
            printf("%d %s creating file %s\n", my_rank, TEST_NAME, filename);
            if ((ret = create_nc_sample_1(iosysid[my_comp_idx], flavor[flv], filename, my_rank, NULL)))
                ERR(ret);

            /* Check the file for correctness. */
            if ((ret = check_nc_sample_1(iosysid[my_comp_idx], flavor[flv], filename, my_rank, NULL)))
                ERR(ret);

        } /* next netcdf format flavor */

        /* If I don't sleep here for a second, there are problems. */
        sleep(2);

        /* Finalize the IO system. Only call this from the computation tasks. */
        if (verbose)
            printf("%d test_intercomm3 Freeing PIO resources\n", my_rank);
        for (int c = 0; c < COMPONENT_COUNT; c++)
        {
            if ((ret = PIOc_finalize(iosysid[c])))
                ERR(ret);
            printf("%d test_intercomm3 PIOc_finalize completed for iosysid = %d\n", my_rank, iosysid[c]);
        }
    } /* endif comp_task */

    /* Wait for everyone to catch up. */
    printf("%d %s waiting for all processes!\n", my_rank, TEST_NAME);
    MPI_Barrier(MPI_COMM_WORLD);

    /* Finalize the MPI library. */
    printf("%d %s Finalizing...\n", my_rank, TEST_NAME);
    if ((ret = pio_test_finalize(&test_comm)))
        return ERR_AWFUL;

    printf("%d %s SUCCESS!!\n", my_rank, TEST_NAME);

    return 0;
}
