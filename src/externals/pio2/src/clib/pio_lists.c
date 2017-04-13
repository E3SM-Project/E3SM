/**
 * @file
 * PIO list functions.
 */
#include <config.h>
#include <pio.h>
#include <pio_internal.h>
#include <string.h>
#include <stdio.h>

static io_desc_t *pio_iodesc_list=NULL;
static io_desc_t *current_iodesc=NULL;
static iosystem_desc_t *pio_iosystem_list=NULL;
static file_desc_t *pio_file_list = NULL;
static file_desc_t *current_file=NULL;

/** Add a new entry to the global list of open files.
 *
 * @param file pointer to the file_desc_t struct for the new file.
 */
void pio_add_to_file_list(file_desc_t *file)
{
    file_desc_t *cfile;

    assert(file);

    /* This file will be at the end of the list, and have no next. */
    file->next = NULL;

    /* Get a pointer to the global list of files. */
    cfile = pio_file_list;

    /* Keep a global pointer to the current file. */
    current_file = file;

    /* If there is nothing in the list, then file will be the first
     * entry. Otherwise, move to end of the list. */
    if (!cfile)
        pio_file_list = file;
    else
    {
        while (cfile->next)
            cfile = cfile->next;
        cfile->next = file;
    }
}

/** Given ncid, find the file_desc_t data for an open file. The ncid
 * used is the interally generated pio_ncid.
 *
 * @param ncid the PIO assigned ncid of the open file.
 * @param cfile1 pointer to a pointer to a file_desc_t. The pointer
 * will get a copy of the pointer to the file info.
 *
 * @returns 0 for success, error code otherwise.
 */
int pio_get_file(int ncid, file_desc_t **cfile1)
{
    file_desc_t *cfile = NULL;

    LOG((2, "pio_get_file ncid = %d", ncid));

    /* Caller must provide this. */
    if (!cfile1)
        return PIO_EINVAL;

    /* Find the file pointer. */
    if (current_file && current_file->pio_ncid == ncid)
        cfile = current_file;
    else
        for (cfile = pio_file_list; cfile; cfile = cfile->next)
            if (cfile->pio_ncid == ncid)
            {
                current_file = cfile;
                break;
            }

    /* If not found, return error. */
    if (!cfile)
        return PIO_EBADID;

    /* We depend on every file having a pointer to the iosystem. */
    if (!cfile->iosystem)
        return PIO_EINVAL;

    /* Let's just ensure we have a valid IO type. */
    pioassert(iotype_is_valid(cfile->iotype), "invalid IO type", __FILE__, __LINE__);

    /* Copy pointer to file info. */
    *cfile1 = cfile;

    return PIO_NOERR;
}

/** Delete a file from the list of open files.
 *
 * @param ncid ID of file to delete from list
 * @returns 0 for success, error code otherwise
 */
int pio_delete_file_from_list(int ncid)
{

    file_desc_t *cfile, *pfile = NULL;

    /* Look through list of open files. */
    for (cfile = pio_file_list; cfile; cfile = cfile->next)
    {
        if (cfile->pio_ncid == ncid)
        {
            if (!pfile)
                pio_file_list = cfile->next;
            else
                pfile->next = cfile->next;

            if (current_file == cfile)
                current_file = pfile;

            /* Free the memory used for this file. */
            free(cfile);
            return PIO_NOERR;
        }
        pfile = cfile;
    }

    /* No file was found. */
    return PIO_EBADID;
}

/** Delete iosystem info from list.
 *
 * @param piosysid the iosysid to delete
 * @returns 0 on success, error code otherwise
 */
int pio_delete_iosystem_from_list(int piosysid)
{
    iosystem_desc_t *ciosystem, *piosystem = NULL;

    LOG((1, "pio_delete_iosystem_from_list piosysid = %d", piosysid));

    for (ciosystem = pio_iosystem_list; ciosystem; ciosystem = ciosystem->next)
    {
        if (ciosystem->iosysid == piosysid)
        {
            if (piosystem == NULL)
                pio_iosystem_list = ciosystem->next;
            else
                piosystem->next = ciosystem->next;
            free(ciosystem);
            return PIO_NOERR;
        }
        piosystem = ciosystem;
    }
    return PIO_EBADID;
}

/** Add iosystem info to list.
 *
 * @param ios pointer to the iosystem_desc_t info to add.
 * @returns 0 on success, error code otherwise
 */
int pio_add_to_iosystem_list(iosystem_desc_t *ios)
{
    iosystem_desc_t *cios;
    int i = 1;

    assert(ios);

    ios->next = NULL;
    cios = pio_iosystem_list;
    if (!cios)
        pio_iosystem_list = ios;
    else
    {
        i++;
        while (cios->next)
        {
            cios = cios->next;
            i++;
        }
        cios->next = ios;
    }

    ios->iosysid = i << 16;

    return ios->iosysid;
}

/** Get iosystem info from list.
 *
 * @param iosysid id of the iosystem
 * @returns pointer to iosystem_desc_t, or NULL if not found.
 */
iosystem_desc_t *pio_get_iosystem_from_id(int iosysid)
{
    iosystem_desc_t *ciosystem;

    LOG((2, "pio_get_iosystem_from_id iosysid = %d", iosysid));

    for (ciosystem = pio_iosystem_list; ciosystem; ciosystem = ciosystem->next)
        if (ciosystem->iosysid == iosysid)
            return ciosystem;

    return NULL;
}

/** Count the number of open iosystems.
 *
 * @param niosysid pointer that will get the number of open iosystems.
 * @returns 0 for success.
 */
int
pio_num_iosystem(int *niosysid)
{
    int count = 0;

    /* Count the elements in the list. */
    for (iosystem_desc_t *c = pio_iosystem_list; c; c = c->next)
        count++;

    /* Return count to caller via pointer. */
    if (niosysid)
        *niosysid = count;

    return PIO_NOERR;
}

/** Add an iodesc.
 *
 * @param io_desc_t pointer to data to add to list.
 * @returns the ioid of the newly added iodesc.
 */
int pio_add_to_iodesc_list(io_desc_t *iodesc)
{
    io_desc_t *ciodesc;
    int imax = 512;

    iodesc->next = NULL;
    if (pio_iodesc_list == NULL)
        pio_iodesc_list = iodesc;
    else{
        imax++;
        for (ciodesc = pio_iodesc_list; ciodesc->next; ciodesc=ciodesc->next, imax=ciodesc->ioid + 1);
        ciodesc->next = iodesc;
    }
    iodesc->ioid = imax;
    current_iodesc = iodesc;

    return iodesc->ioid;
}

/** Get an iodesc.
 *
 * @param ioid ID of iodesc to get.
 * @returns pointer to the iodesc struc.
 */
io_desc_t *pio_get_iodesc_from_id(int ioid)
{
    io_desc_t *ciodesc = NULL;

    if (current_iodesc != NULL && current_iodesc->ioid == abs(ioid))
        ciodesc = current_iodesc;
    else
        for (ciodesc = pio_iodesc_list; ciodesc; ciodesc = ciodesc->next)
            if (ciodesc->ioid == abs(ioid))
            {
                current_iodesc = ciodesc;
                break;
            }

    return ciodesc;
}

/** Delete an iodesc.
 *
 * @param ioid ID of iodesc to delete.
 * @returns 0 on success, error code otherwise.
 */
int pio_delete_iodesc_from_list(int ioid)
{
    io_desc_t *ciodesc, *piodesc = NULL;

    for (ciodesc = pio_iodesc_list; ciodesc; ciodesc = ciodesc->next)
    {
        if (ciodesc->ioid == ioid)
        {
            if (piodesc == NULL)
                pio_iodesc_list = ciodesc->next;
            else
                piodesc->next = ciodesc->next;

            if (current_iodesc == ciodesc)
                current_iodesc = pio_iodesc_list;
            free(ciodesc);
            return PIO_NOERR;
        }
        piodesc = ciodesc;
    }
    return PIO_EBADID;
}
