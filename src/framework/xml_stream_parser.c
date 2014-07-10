#include <stdio.h>
#include <stdlib.h>
#include "ezxml/ezxml.h"
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>


/* 
 * Interface routines for building streams at run-time; defined in mpas_stream_manager.F
 */
void stream_mgr_create_stream_c(void *, const char *, int *);


int par_read(char *fname, int *mpi_comm, char **xml_buf, size_t *bufsize)
{
	int iofd;
	int rank;
	struct stat s;
	int err;

#ifdef _MPI
#include "mpi.h"
	
	MPI_Comm comm;

	comm = MPI_Comm_f2c((MPI_Fint)(*mpi_comm));
	err = MPI_Comm_rank(comm, &rank);
#else
	rank = 0;
#endif

	if (rank == 0) {
		iofd = open(fname, O_RDONLY);
		if (!iofd) {
			fprintf(stderr, "Error: Could not open run-time I/O config file %s\n", fname);
			return 1;
		}

		fstat(iofd, &s);
		*bufsize = (size_t)s.st_size;
#ifdef _MPI
		err = MPI_Bcast((void *)bufsize, (int)sizeof(size_t), MPI_BYTE, 0, comm);
#endif
	
		*xml_buf = (char *)malloc(*bufsize);
		read(iofd, (void *)(*xml_buf), *bufsize);

#ifdef _MPI
		err = MPI_Bcast((void *)(*xml_buf), (int)(*bufsize), MPI_CHAR, 0, comm);
#endif
	}
	else {
#ifdef _MPI
		err = MPI_Bcast((void *)bufsize, (int)sizeof(size_t), MPI_BYTE, 0, comm);
#endif
		*xml_buf = (char *)malloc(*bufsize);

#ifdef _MPI
		err = MPI_Bcast((void *)(*xml_buf), (int)(*bufsize), MPI_CHAR, 0, comm);
#endif
	}

	return 0;
}


void xml_stream_parser(char *fname, void *manager, int *mpi_comm)
{
        char *xml_buf;
	size_t bufsize;
	ezxml_t streams;
	ezxml_t foo_xml;
	const char *attr;
	int err;

	fprintf(stderr, "MGD DEV begin parsing run-time I/O from %s\n", fname);

        if (par_read(fname, mpi_comm, &xml_buf, &bufsize) != 0) {
		return;
	}

	streams = ezxml_parse_str(xml_buf, bufsize);
	if (!streams) {
		fprintf(stderr, "Error: Problems encountered while parsing run-time I/O config file %s\n", fname);
		return;
	}	

	err = 0;
	for (foo_xml = ezxml_child(streams, "stream"); foo_xml; foo_xml = ezxml_next(foo_xml)) {
		attr = ezxml_attr(foo_xml, "name");
		fprintf(stderr, "MGD DEV Found stream named %s\n", attr);
		stream_mgr_create_stream_c(manager, attr, &err);
		err++;
	}

	free(xml_buf);

	fprintf(stderr, "MGD DEV done parsing run-time I/O\n");
}
