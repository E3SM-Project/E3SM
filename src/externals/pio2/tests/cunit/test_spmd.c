/*
 * This program tests some internal functions in the library.
 *
 * Jim Edwards
 * Ed Hartnett, 11/23/16
 */
#include <pio.h>
#include <pio_tests.h>
#include <pio_internal.h>

/* The number of tasks this test should run on. */
#define TARGET_NTASKS 4

/* The minimum number of tasks this test should run on. */
#define MIN_NTASKS 1

/* The name of this test. */
#define TEST_NAME "test_spmd"

/* Number of test cases in inner loop of test. */
#define NUM_TEST_CASES 5

#define TEST_MAX_GATHER_BLOCK_SIZE 32

/* Test MPI_Alltoallw by having processor i send different amounts of
 * data to each processor.  The first test sends i items to processor
 * i from all processors. */
int run_spmd_tests(MPI_Comm test_comm)
{
    int my_rank;  /* 0-based rank in test_comm. */
    int ntasks;   /* Number of tasks in test_comm. */
    int num_elem; /* Number of elements in buffers. */
    int type_size; /* Size in bytes of an element. */
    int mpierr;   /* Return value from MPI calls. */
    int ret;      /* Return value. */

    /* Learn rank and size. */
    if ((mpierr = MPI_Comm_size(test_comm, &ntasks)))
        MPIERR(mpierr);
    if ((mpierr = MPI_Comm_rank(test_comm, &my_rank)))
        MPIERR(mpierr);

    /* Determine size of buffers. */
    num_elem = ntasks;

    int sbuf[ntasks];    /* The send buffer. */
    int rbuf[ntasks];    /* The receive buffer. */
    int sendcounts[ntasks]; /* Number of elements of data being sent from each task. */
    int recvcounts[ntasks]; /* Number of elements of data being sent from each task. */
    int sdispls[ntasks]; /* Displacements for sending data. */
    int rdispls[ntasks]; /* Displacements for receiving data. */
    MPI_Datatype sendtypes[ntasks]; /* MPI types of data being sent. */
    MPI_Datatype recvtypes[ntasks]; /* MPI types of data being received. */

    /* Load up the send buffer. */
    for (int i = 0; i < num_elem; i++)
        sbuf[i] = my_rank;

    /* Load up the receive buffer to make debugging easier. */
    for (int i = 0; i < num_elem; i++)
        rbuf[i] = -999;

    /* Get the size of the int type for MPI. (Should always be 4.) */
    if ((mpierr = MPI_Type_size(MPI_INT, &type_size)))
        return check_mpi(NULL, mpierr, __FILE__, __LINE__);
    assert(type_size == sizeof(int));

    /* Initialize the arrays. */
    for (int i = 0; i < ntasks; i++)
    {
        sendcounts[i] = 1;
        sdispls[i] = 0;
        sendtypes[i] = MPI_INT;
        recvcounts[i] = 1;
        rdispls[i] = i * type_size;
        recvtypes[i] = MPI_INT;
    }

    /* Perform tests for different values of msg_cnt. (BTW it hangs
     * with msg_cnt = 1!). */
    for (int msg_cnt = 0; msg_cnt < TARGET_NTASKS; msg_cnt = msg_cnt ? msg_cnt * 2 : 4)
    {
        if (!my_rank)
            printf("message count %d\n",msg_cnt);

        for (int itest = 0; itest < NUM_TEST_CASES; itest++)
        {
            bool hs = false;
            bool isend = false;

            /* Wait for all tasks. */
            MPI_Barrier(test_comm);

            /* Print results. */
            if (!my_rank)
                for (int e = 0; e < num_elem; e++)
                    printf("sbuf[%d] = %d\n", e, sbuf[e]);

            /* Set the parameters different for each test case. */
            if (itest == 1)
            {
                hs = true;
                isend = true;
            }
            else if (itest == 2)
            {
                hs = false;
                isend = true;
            }
            else if (itest == 3)
            {
                hs = false;
                isend = false;
            }
            else if (itest == 4)
            {
                hs = true;
                isend = false;
            }

            /* Run the swapm function. */
            if ((ret = pio_swapm(sbuf, sendcounts, sdispls, sendtypes, rbuf, recvcounts,
                                 rdispls, recvtypes, test_comm, hs, isend, msg_cnt)))
                return ret;

            /* Print results. */
            /* MPI_Barrier(test_comm); */
            /* for (int e = 0; e < num_elem; e++) */
            /*     printf("%d sbuf[%d] = %d\n", my_rank, e, sbuf[e]); */
            /* MPI_Barrier(test_comm); */
            /* for (int e = 0; e < num_elem; e++) */
            /*     printf("%d rbuf[%d] = %d\n", my_rank, e, rbuf[e]); */

            /* Check that rbuf has 0, 1, ..., ntasks-1. */
            for (int e = 0; e < num_elem; e++)
                if (((int *)rbuf)[e] != e)
                    return ERR_WRONG;
        }
    }

    return 0;
}

/* Test some of the functions in the file pioc_sc.c. 
 *
 * @param test_comm the MPI communicator that the test code is running on. 
 * @returns 0 for success, error code otherwise.
 */
int run_sc_tests(MPI_Comm test_comm)
{
#define SC_ARRAY_LEN 3
    int my_rank;  /* 0-based rank in test_comm. */
    int ntasks;   /* Number of tasks in test_comm. */
    int mpierr;   /* Return value from MPI calls. */
    int array1[SC_ARRAY_LEN] = {7, 42, 14};
    int array2[SC_ARRAY_LEN] = {2, 3, 7};
    int array3[SC_ARRAY_LEN] = {90, 180, 270};
    int array4[SC_ARRAY_LEN] = {1, 180, 270};

    /* Learn rank and size. */
    if ((mpierr = MPI_Comm_size(test_comm, &ntasks)))
        MPIERR(mpierr);
    if ((mpierr = MPI_Comm_rank(test_comm, &my_rank)))
        MPIERR(mpierr);

    /* Test the gcd() function. */
    if (gcd(0, 2) != 2)
        return ERR_WRONG;
    if (gcd(2, 2) != 2)
        return ERR_WRONG;
    if (gcd(42, 2) != 2)
        return ERR_WRONG;

    /* Test the long long version. */
    if (lgcd(0, 2) != 2)
        return ERR_WRONG;
    if (lgcd(2, 2) != 2)
        return ERR_WRONG;
    if (lgcd(42, 2) != 2)
        return ERR_WRONG;

    /* Test the gcd_array() function. */
    if (gcd_array(SC_ARRAY_LEN, array1) != 7)
        return ERR_WRONG;
    if (gcd_array(SC_ARRAY_LEN, array2) != 1)
        return ERR_WRONG;
    if (gcd_array(SC_ARRAY_LEN, array3) != 90)
        return ERR_WRONG;
    if (gcd_array(SC_ARRAY_LEN, array4) != 1)
        return ERR_WRONG;

    return 0;
}

/* This test code was recovered from main() in pioc_sc.c. */
int test_CalcStartandCount()
{
    int ndims = 2;
    int gdims[2] = {31, 777602};
    int num_io_procs = 24;
    bool converged = false;
    PIO_Offset start[ndims], kount[ndims];
    int iorank, numaiotasks = 0;
    long int tpsize = 0;
    long int psize;
    long int pgdims = 1;
    int scnt;

    for (int i = 0; i < ndims; i++)
        pgdims *= gdims[i];

    while (!converged)
    {
        for (iorank = 0; iorank < num_io_procs; iorank++)
        {
            numaiotasks = CalcStartandCount(PIO_DOUBLE, ndims, gdims, num_io_procs, iorank,
                                            start, kount);
            if (iorank < numaiotasks)
                printf("iorank %d start %lld %lld count %lld %lld\n", iorank, start[0],
                       start[1], kount[0], kount[1]);

            if (numaiotasks < 0)
                return numaiotasks;

            psize = 1;
            scnt = 0;
            for (int i = 0; i < ndims; i++)
            {
                psize *= kount[i];
                scnt += kount[i];
            }
            tpsize += psize;
        }

        if (tpsize == pgdims)
            converged = true;
        else
        {
            printf("Failed to converge %ld %ld %d\n", tpsize, pgdims, num_io_procs);
            tpsize = 0;
            num_io_procs--;
        }
    }

    return 0;
}

/* Tesst some list stuff. */
int test_lists()
{
    file_desc_t *fdesc;
    
    /* Test that bad input is correctly rejected. */
    if (pio_delete_iodesc_from_list(42) != PIO_EBADID)
        return ERR_WRONG;
    if (pio_delete_iosystem_from_list(42) != PIO_EBADID)
        return ERR_WRONG;
    if (pio_delete_file_from_list(42) != PIO_EBADID)
        return ERR_WRONG;
    if (pio_get_file(42, NULL) != PIO_EINVAL)
        return ERR_WRONG;
    if (pio_get_file(42, &fdesc) != PIO_EBADID)
        return ERR_WRONG;
    return 0;
}

/* Test some of the rearranger utility functions. */
int test_rearranger_opts1()
{
    rearr_comm_fc_opt_t *ro1;
    rearr_comm_fc_opt_t *ro2;
    rearr_comm_fc_opt_t *ro3;

    if (!(ro1 = calloc(1, sizeof(rearr_comm_fc_opt_t))))
        return ERR_AWFUL;
    if (!(ro2 = calloc(1, sizeof(rearr_comm_fc_opt_t))))
        return ERR_AWFUL;
    if (!(ro3 = calloc(1, sizeof(rearr_comm_fc_opt_t))))
        return ERR_AWFUL;

    /* This should not work. */
    if (PIOc_set_rearr_opts(42, 1, 1, 0, 0, 0, 0, 0, 0) != PIO_EBADID)
        return ERR_WRONG;

    /* ro1 and ro2 are the same. */
    if (!cmp_rearr_comm_fc_opts(ro1, ro2))
        return ERR_WRONG;

    /* Make ro3 different. */
    ro3->enable_hs = 1;
    if (cmp_rearr_comm_fc_opts(ro1, ro3))
        return ERR_WRONG;
    ro3->enable_hs = 0;
    ro3->enable_isend = 1;
    if (cmp_rearr_comm_fc_opts(ro1, ro3))
        return ERR_WRONG;
    ro3->enable_isend = 0;
    ro3->max_pend_req = 1;
    if (cmp_rearr_comm_fc_opts(ro1, ro3))
        return ERR_WRONG;

    /* Free resourses. */
    free(ro1);
    free(ro2);
    free(ro3);
    
    return 0;
}

/* Test some of the rearranger utility functions. */
int test_rearranger_opts2()
{
    iosystem_desc_t my_ios;
    iosystem_desc_t *ios = &my_ios;

    /* I'm not sure what the point of this function is... */
    check_and_reset_rearr_opts(ios);
    
    return 0;
}

/* Test the compare_offsets() function. */
int test_compare_offsets()
{
    mapsort m1, m2, m3;

    m1.rfrom = 0;
    m1.soffset = 0;
    m1.iomap = 0;
    m2.rfrom = 0;
    m2.soffset = 0;
    m2.iomap = 0;
    m3.rfrom = 0;
    m3.soffset = 0;
    m3.iomap = 1;

    /* Return 0 if either or both parameters are null. */
    if (compare_offsets(NULL, &m2))
        return ERR_WRONG;
    if (compare_offsets(&m1, NULL))
        return ERR_WRONG;
    if (compare_offsets(NULL, NULL))
        return ERR_WRONG;

    /* m1 and m2 are the same. */
    if (compare_offsets(&m1, &m2))
        return ERR_WRONG;

    /* m1 and m3 are different. */
    if (compare_offsets(&m1, &m3) != -1)
        return ERR_WRONG;
    return 0;
}

/* Test the ceil2() and pair() functions. */
int test_ceil2_pair()
{
    /* Test the ceil2() function. */
    if (ceil2(1) != 1)
        return ERR_WRONG;
    if (ceil2(-100) != 1)
        return ERR_WRONG;
    if (ceil2(2) != 2)
        return ERR_WRONG;
    if (ceil2(3) != 4)
        return ERR_WRONG;
    if (ceil2(16) != 16)
        return ERR_WRONG;
    if (ceil2(17) != 32)
        return ERR_WRONG;

    /* Test the pair() function. */
    if (pair(4, 0, 0) != 1)
        return ERR_WRONG;
    if (pair(4, 2, 2) != 1)
        return ERR_WRONG;
    
    return 0;
}

/* Test the function that finds an MPI type to match a PIO type. */
int test_find_mpi_type()
{
    MPI_Datatype mpi_type;
    int ret;

    /* This should not work. */
    if (find_mpi_type(PIO_BYTE + 42, &mpi_type) != PIO_EBADTYPE)
        return ERR_WRONG;

    /* Try every atomic type. */
    if ((ret = find_mpi_type(PIO_BYTE, &mpi_type)))
        return ret;
    if (mpi_type != MPI_BYTE)
        return ERR_WRONG;

    if ((ret = find_mpi_type(PIO_CHAR, &mpi_type)))
        return ret;
    if (mpi_type != MPI_CHAR)
        return ERR_WRONG;

    if ((ret = find_mpi_type(PIO_SHORT, &mpi_type)))
        return ret;
    if (mpi_type != MPI_SHORT)
        return ERR_WRONG;

    if ((ret = find_mpi_type(PIO_INT, &mpi_type)))
        return ret;
    if (mpi_type != MPI_INT)
        return ERR_WRONG;

    if ((ret = find_mpi_type(PIO_FLOAT, &mpi_type)))
        return ret;
    if (mpi_type != MPI_FLOAT)
        return ERR_WRONG;

    if ((ret = find_mpi_type(PIO_DOUBLE, &mpi_type)))
        return ret;
    if (mpi_type != MPI_DOUBLE)
        return ERR_WRONG;

#ifdef _NETCDF4
    if ((ret = find_mpi_type(PIO_UBYTE, &mpi_type)))
        return ret;
    if (mpi_type != MPI_UNSIGNED_CHAR)
        return ERR_WRONG;

    if ((ret = find_mpi_type(PIO_USHORT, &mpi_type)))
        return ret;
    if (mpi_type != MPI_UNSIGNED_SHORT)
        return ERR_WRONG;

    if ((ret = find_mpi_type(PIO_UINT, &mpi_type)))
        return ret;
    if (mpi_type != MPI_UNSIGNED)
        return ERR_WRONG;

    if ((ret = find_mpi_type(PIO_INT64, &mpi_type)))
        return ret;
    if (mpi_type != MPI_LONG_LONG)
        return ERR_WRONG;

    if ((ret = find_mpi_type(PIO_UINT64, &mpi_type)))
        return ret;
    if (mpi_type != MPI_UNSIGNED_LONG_LONG)
        return ERR_WRONG;

    if ((ret = find_mpi_type(PIO_STRING, &mpi_type)))
        return ret;
    if (mpi_type != MPI_CHAR)
        return ERR_WRONG;

#endif /* _NETCDF4 */
    return PIO_NOERR;
}

int test_misc()
{
    wmulti_buffer wmb;

    /* This should not work. */
    if (flush_buffer(TEST_VAL_42, &wmb, 0) != PIO_EBADID)
        return ERR_WRONG;
    
    return 0;
}

/* Run Tests for pio_spmd.c functions. */
int main(int argc, char **argv)
{
    int my_rank; /* Zero-based rank of processor. */
    int ntasks;  /* Number of processors involved in current execution. */
    int ret;     /* Return code. */
    MPI_Comm test_comm; /* A communicator for this test. */

    /* Initialize test. */
    if ((ret = pio_test_init2(argc, argv, &my_rank, &ntasks, MIN_NTASKS,
                              TARGET_NTASKS, 3, &test_comm)))
        ERR(ERR_INIT);

    /* Test code runs on TARGET_NTASKS tasks. The left over tasks do
     * nothing. */
    if (my_rank < TARGET_NTASKS)
    {
        printf("%d running tests for functions in pioc_sc.c\n", my_rank);
        if ((ret = run_sc_tests(test_comm)))
            return ret;

        printf("%d running spmd test code\n", my_rank);
        if ((ret = run_spmd_tests(test_comm)))
            return ret;
        
        printf("%d running CalcStartandCount test code\n", my_rank);
        if ((ret = test_CalcStartandCount()))
            return ret;

        printf("%d running list tests\n", my_rank);
        if ((ret = test_lists()))
            return ret;

        printf("%d running rearranger opts tests 1\n", my_rank);
        if ((ret = test_rearranger_opts1()))
            return ret;

        printf("%d running rearranger opts tests 2\n", my_rank);
        if ((ret = test_rearranger_opts2()))
            return ret;

        printf("%d running compare_offsets tests\n", my_rank);
        if ((ret = test_compare_offsets()))
            return ret;

        printf("%d running ceil2/pair tests\n", my_rank);
        if ((ret = test_ceil2_pair()))
            return ret;

        printf("%d running find_mpi_type tests\n", my_rank);
        if ((ret = test_find_mpi_type()))
            return ret;

        printf("%d running misc tests\n", my_rank);
        if ((ret = test_misc()))
            return ret;

    } /* endif my_rank < TARGET_NTASKS */

    /* Finalize the MPI library. */
    printf("%d %s Finalizing...\n", my_rank, TEST_NAME);
    if ((ret = pio_test_finalize(&test_comm)))
        return ret;

    printf("%d %s SUCCESS!!\n", my_rank, TEST_NAME);

    return 0;
}
