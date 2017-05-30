/*
 * Tests the PIO library with multiple iosysids in use at the
 * same time.
 *
 * This is a simplified, C version of the fortran
 * pio_iosystem_tests3.F90.
 *
 * Ed Hartnett
 */
#include <pio.h>
#include <pio_tests.h>

/* The number of tasks this test should run on. */
#define TARGET_NTASKS 4

/* The name of this test. */
#define TEST_NAME "test_iosystem3"

/* Used to define netcdf test file. */
#define PIO_TF_MAX_STR_LEN 100
#define ATTNAME "filename"
#define DIMNAME "filename_dim"

/* Used to devide up the tasks into MPI groups. */
#define OVERLAP_NUM_RANGES 2
#define EVEN_NUM_RANGES 1

/* Used when initializing PIO. */
#define STRIDE1 1
#define STRIDE2 2
#define BASE0 0
#define BASE1 1
#define NUM_IO1 1
#define NUM_IO2 2
#define NUM_IO4 4
#define REARRANGER 2

/* This creates a netCDF file in the specified format, with some
 * sample values. */
int create_file(MPI_Comm comm, int iosysid, int format, char *filename,
                char *attname, char *dimname, int my_rank)
{
    int ncid, varid, dimid;
    int ret;

    /* Create the file. */
    if ((ret = PIOc_createfile(iosysid, &ncid, &format, filename, NC_CLOBBER)))
        return ret;
    printf("%d file created ncid = %d\n", my_rank, ncid);

    /* Define a dimension. */
    printf("%d defining dimension %s\n", my_rank, dimname);
    if ((ret = PIOc_def_dim(ncid, dimname, PIO_TF_MAX_STR_LEN, &dimid)))
        return ret;

    /* Define a 1-D variable. */
    printf("%d defining variable %s\n", my_rank, attname);
    if ((ret = PIOc_def_var(ncid, attname, NC_CHAR, 1, &dimid, &varid)))
        return ret;

    /* Write an attribute. */
    if ((ret = PIOc_put_att_text(ncid, varid, attname, strlen(filename), filename)))
        return ret;

    /* End define mode. */
    printf("%d ending define mode ncid = %d\n", my_rank, ncid);
    if ((ret = PIOc_enddef(ncid)))
        return ret;
    printf("%d define mode ended ncid = %d\n", my_rank, ncid);

    /* Close the file. */
    printf("%d closing file ncid = %d\n", my_rank, ncid);
    if ((ret = PIOc_closefile(ncid)))
        return ret;
    printf("%d closed file ncid = %d\n", my_rank, ncid);

    return PIO_NOERR;
}

/* This checks an already-open netCDF file. */
int check_file(MPI_Comm comm, int iosysid, int format, int ncid, char *filename,
               char *attname, char *dimname, int my_rank)
{
    int dimid;
    int varid;
    char *att_data;
    int ret;

    /* Check the dimid. */
    if ((ret = PIOc_inq_dimid(ncid, dimname, &dimid)))
        return ret;
    if (dimid)
        return ERR_WRONG;

    /* Check the varid (it's got the same name as the att). */
    if ((ret = PIOc_inq_varid(ncid, attname, &varid)))
        return ret;
    if (varid)
        return ERR_WRONG;

    /* Check the attribute. Null terminating byte deliberately ignored
     * to match fortran code. */
    if (!(att_data = malloc(strlen(filename) * sizeof(char))))
        return PIO_ENOMEM;
    if ((ret = PIOc_get_att(ncid, varid, attname, att_data)))
        return ret;
    printf("%d DONE with get_att!!!\n", my_rank);
    if (strncmp(att_data, filename, strlen(filename)))
        return ERR_WRONG;
    free(att_data);
    printf("%d DONE with get_att!!!\n", my_rank);

    return PIO_NOERR;
}

/* This opens and checks a netCDF file. */
int open_and_check_file(MPI_Comm comm, int iosysid, int iotype, int *ncid, char *fname,
                        char *attname, char *dimname, int disable_close, int my_rank)
{
    int mode = PIO_WRITE;
    int ret;

    /* Open the file. */
    if ((ret = PIOc_openfile(iosysid, ncid, &iotype, fname, mode)))
        return ret;

    /* Check the file. */
    if ((ret = check_file(comm, iosysid, iotype, *ncid, fname, attname, dimname, my_rank)))
        return ret;

    /* Close the file, maybe. */
    if (!disable_close)
        if ((ret = PIOc_closefile(*ncid)))
            return ret;

    return PIO_NOERR;
}

/* Run async tests. */
int main(int argc, char **argv)
{
    int my_rank; /* Zero-based rank of processor. */
    int ntasks; /* Number of processors involved in current execution. */
    int iosysid_world; /* The ID for the parallel I/O system. */
    int even_iosysid; /* The ID for iosystem of even_comm. */
    int overlap_iosysid; /* The ID for iosystem of even_comm. */
    MPI_Group world_group; /* An MPI group of world. */
    MPI_Group even_group; /* An MPI group of 0 and 2. */
    MPI_Group overlap_group; /* An MPI group of 0, 1, and 3. */
    MPI_Comm even_comm = MPI_COMM_NULL; /* Communicator for tasks 0, 2 */
    MPI_Comm overlap_comm = MPI_COMM_NULL; /* Communicator for tasks 0, 1, 2. */
    int even_rank = -1, overlap_rank = -1; /* Tasks rank in communicator. */
    int even_size = 0, overlap_size = 0; /* Size of communicator. */
    int num_flavors; /* Number of PIO netCDF flavors in this build. */
    int flavor[NUM_FLAVORS]; /* iotypes for the supported netCDF IO flavors. */
    int ret; /* Return code. */
    MPI_Comm test_comm;

    /* Initialize test. */
    if ((ret = pio_test_init(argc, argv, &my_rank, &ntasks, TARGET_NTASKS,
                             &test_comm)))
        ERR(ERR_INIT);

    /* Test code runs on TARGET_NTASKS tasks. The left over tasks do
     * nothing. */
    if (my_rank < TARGET_NTASKS)
    {
        /* Figure out iotypes. */
        if ((ret = get_iotypes(&num_flavors, flavor)))
            ERR(ret);

        /* Initialize PIO system on world. */
        printf("%d about to call Init_Intracomm\n", my_rank);
        if ((ret = PIOc_Init_Intracomm(test_comm, NUM_IO4, STRIDE1, BASE0, REARRANGER, &iosysid_world)))
            ERR(ret);
        printf("%d done with Init_Intracomm\n", my_rank);

        /* Set the error handler. */
        /*PIOc_Set_IOSystem_Error_Handling(iosysid_world, PIO_BCAST_ERROR);*/
        printf("%d about to set iosystem error hanlder for world\n", my_rank);
        if ((ret = PIOc_set_iosystem_error_handling(iosysid_world, PIO_BCAST_ERROR, NULL)))
            ERR(ret);
        printf("%d done setting iosystem error hanlder for world\n", my_rank);

        /* Get MPI_Group of world comm. */
        if ((ret = MPI_Comm_group(test_comm, &world_group)))
            ERR(ret);

        /* Create a group with tasks 0 and 2. */
        int even_ranges[EVEN_NUM_RANGES][3] = {{0, 2, 2}};
        if ((ret = MPI_Group_range_incl(world_group, EVEN_NUM_RANGES, even_ranges, &even_group)))
            ERR(ret);

        /* Create a communicator from the even_group. */
        if ((ret = MPI_Comm_create(test_comm, even_group, &even_comm)))
            ERR(ret);

        /* Learn my rank and the total number of processors in even group. */
        if (even_comm != MPI_COMM_NULL)
        {
            if ((ret = MPI_Comm_rank(even_comm, &even_rank)))
                MPIERR(ret);
            if ((ret = MPI_Comm_size(even_comm, &even_size)))
                MPIERR(ret);
        }
        printf("%d even_comm = %d even_rank = %d even_size = %d\n", my_rank,
               even_comm, even_rank, even_size);

        /* Create a group with tasks 0, 1, and 3. */
        int overlap_ranges[OVERLAP_NUM_RANGES][3] = {{0, 0, 1}, {1, 3, 2}};
        if ((ret = MPI_Group_range_incl(world_group, OVERLAP_NUM_RANGES,
                                        overlap_ranges, &overlap_group)))
            ERR(ret);

        /* Create a communicator from the overlap_group. */
        if ((ret = MPI_Comm_create(test_comm, overlap_group, &overlap_comm)))
            ERR(ret);

        /* Learn my rank and the total number of processors in overlap
         * group. */
        if (overlap_comm != MPI_COMM_NULL)
        {
            if ((ret = MPI_Comm_rank(overlap_comm, &overlap_rank)))
                MPIERR(ret);
            if ((ret = MPI_Comm_size(overlap_comm, &overlap_size)))
                MPIERR(ret);
        }
        printf("%d overlap_comm = %d overlap_rank = %d overlap_size = %d\n", my_rank,
               overlap_comm, overlap_rank, overlap_size);

        /* Initialize PIO system for even. */
        if (even_comm != MPI_COMM_NULL)
        {
            if ((ret = PIOc_Init_Intracomm(even_comm, NUM_IO1, STRIDE1, BASE1, REARRANGER, &even_iosysid)))
                ERR(ret);

            /* These should not work. */
            if (PIOc_set_hint(even_iosysid + TEST_VAL_42, NULL, NULL) != PIO_EBADID)
                ERR(ERR_WRONG);
            if (PIOc_set_hint(even_iosysid, NULL, NULL) != PIO_EINVAL)
                ERR(ERR_WRONG);

            /* Set the hint (which will be ignored). */
            if ((ret = PIOc_set_hint(even_iosysid, "hint", "hint_value")))
                ERR(ret);

            /* Set the error handler. */
            /*PIOc_Set_IOSystem_Error_Handling(even_iosysid, PIO_BCAST_ERROR);*/
            printf("%d about to set iosystem error hanlder for even\n", my_rank);
            if ((ret = PIOc_set_iosystem_error_handling(even_iosysid, PIO_BCAST_ERROR, NULL)))
                ERR(ret);
            printf("%d done setting iosystem error hanlder for even\n", my_rank);
        }

        /* Initialize PIO system for overlap comm. */
        if (overlap_comm != MPI_COMM_NULL)
        {
            if ((ret = PIOc_Init_Intracomm(overlap_comm, NUM_IO2, STRIDE1, BASE1, REARRANGER,
                                           &overlap_iosysid)))
                ERR(ret);

            printf("%d about to set iosystem error hanlder for overlap\n", my_rank);
            /* Set the error handler. */
            /* if ((ret = PIOc_set_iosystem_error_handling(overlap_iosysid, PIO_BCAST_ERROR))) */
            /*     ERR(ret); */
            PIOc_Set_IOSystem_Error_Handling(overlap_iosysid, PIO_BCAST_ERROR);
            printf("%d done setting iosystem error hanlder for overlap\n", my_rank);
        }

        for (int i = 0; i < num_flavors; i++)
        {
            char fname0[] = "pio_iosys_test_file0.nc";
            char fname1[] = "pio_iosys_test_file1.nc";
            char fname2[] = "pio_iosys_test_file2.nc";
            printf("\n\n%d i = %d\n", my_rank, i);

            if ((ret = create_file(test_comm, iosysid_world, flavor[i], fname0, ATTNAME,
                                   DIMNAME, my_rank)))
                ERR(ret);

            if ((ret = create_file(test_comm, iosysid_world, flavor[i], fname1, ATTNAME,
                                   DIMNAME, my_rank)))
                ERR(ret);

            if ((ret = create_file(test_comm, iosysid_world, flavor[i], fname2, ATTNAME,
                                   DIMNAME, my_rank)))
                ERR(ret);

            /* Now check the first file from WORLD communicator. */
            int ncid;
            if ((ret = open_and_check_file(test_comm, iosysid_world, flavor[i], &ncid, fname0,
                                           ATTNAME, DIMNAME, 1, my_rank)))
                ERR(ret);

            /* Now have the even communicators check the files. */
            int ncid2;
            if (even_comm != MPI_COMM_NULL)
            {
                printf("\n***\n%d Checking file for even_comm\n", my_rank);
                if ((ret = open_and_check_file(even_comm, even_iosysid, flavor[i], &ncid2, fname2,
                                               ATTNAME, DIMNAME, 1, my_rank)))
                    ERR(ret);
                if ((ret = check_file(even_comm, even_iosysid, flavor[i], ncid2, fname2,
                                      ATTNAME, DIMNAME, my_rank)))
                    ERR(ret);
            }

            /* Now have the overlap communicators check the files. */
            int ncid3;
            if (overlap_comm != MPI_COMM_NULL)
            {
                printf("\n***%d Checking file for overlap_comm\n", my_rank);
                if ((ret = open_and_check_file(overlap_comm, overlap_iosysid, flavor[i], &ncid3, fname1,
                                               ATTNAME, DIMNAME, 1, my_rank)))
                    ERR(ret);
                if ((ret = check_file(overlap_comm, overlap_iosysid, flavor[i], ncid3, fname1,
                                      ATTNAME, DIMNAME, my_rank)))
                    ERR(ret);
            }

            /* Close the still-open files. */
            if ((ret = PIOc_closefile(ncid)))
                ERR(ret);
            if (even_comm != MPI_COMM_NULL)
            {
                if ((ret = PIOc_closefile(ncid2)))
                    ERR(ret);
            }
            if (overlap_comm != MPI_COMM_NULL)
            {
                if ((ret = PIOc_closefile(ncid3)))
                    ERR(ret);
            }
        } /* next iotype */
        /* Finalize PIO systems. */
        printf("%d pio finalizing %d\n", my_rank, even_iosysid);
        if (even_comm != MPI_COMM_NULL)
            if ((ret = PIOc_finalize(even_iosysid)))
                ERR(ret);
        printf("%d pio finalizing %d\n", my_rank, overlap_iosysid);
        if (overlap_comm != MPI_COMM_NULL)
        {
            printf("%d calling PIOc_finalize with iosysid = %d\n", my_rank, overlap_iosysid);
            if ((ret = PIOc_finalize(overlap_iosysid)))
                ERR(ret);
        }
        printf("%d pio finalized\n", my_rank);
        if ((ret = PIOc_finalize(iosysid_world)))
            ERR(ret);

        /* Free MPI resources used by test. */
        if ((ret = MPI_Group_free(&overlap_group)))
            ERR(ret);
        if ((ret = MPI_Group_free(&even_group)))
            ERR(ret);
        if ((ret = MPI_Group_free(&world_group)))
            ERR(ret);
        if (overlap_comm != MPI_COMM_NULL)
            if ((ret = MPI_Comm_free(&overlap_comm)))
                ERR(ret);
        if (even_comm != MPI_COMM_NULL)
            if ((ret = MPI_Comm_free(&even_comm)))
                ERR(ret);

    } /* my_rank < TARGET_NTASKS */

    /* Finalize test. */
    printf("%d %s finalizing...\n", my_rank, TEST_NAME);
    if ((ret = pio_test_finalize(&test_comm)))
        return ERR_AWFUL;

    printf("%d %s SUCCESS!!\n", my_rank, TEST_NAME);

    return 0;
}
