// Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
// and the University Corporation for Atmospheric Research (UCAR).
//
// Unless noted otherwise source code is licensed under the BSD license.
// Additional copyright and license information can be found in the LICENSE file
// distributed with this code, or at http://mpas-dev.github.com/license.html
//

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


/*********************************************************************************
 *
 *  Function: par_read
 *
 *  Reads the contents of a file into a buffer in distributed-memory parallel code.
 * 
 *  The buffer xml_buf is allocated with size bufsize, which will be exactly the 
 *  number of bytes in the file fname. Only the master task will actually read the
 *  file, and the contents are broadcast to all other tasks. The mpi_comm argument
 *  is a Fortran MPI communicator used to determine which task is the master task.
 * 
 *  A return code of 0 indicates the file was successfully read and broadcast to
 *  all MPI tasks that belong to the communicator.
 *
 *********************************************************************************/
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


/*********************************************************************************
 *
 *  Function: xml_stream_parser
 *
 *  Parses an XML file and builds streams using the MPAS_stream_manager module
 *  based on the contents of the file.
 *
 *  The fname argument provides the name of the XML file that contains the stream
 *  definitions, manager is a Fortran derived type used by the stream mananger,
 *  and mpi_comm is the Fortran MPI communicator used by MPAS.
 *
 *********************************************************************************/
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
