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
#include <string.h>


/* 
 * Interface routines for building streams at run-time; defined in mpas_stream_manager.F
 */
void stream_mgr_create_stream_c(void *, const char *, int *);


/* 
 * Global variables
 */
static char *global_file;

#define MSGSIZE 128


/*********************************************************************************
 *
 *  Function: fmt_err
 *
 *  Prints an error message in a standard format.
 *
 *********************************************************************************/
void fmt_err(const char *mesg)
{
	fprintf(stderr,"********************************************************************************\n");
	fprintf(stderr,"Error: In file %s, %s\n", global_file, mesg);
	fprintf(stderr,"********************************************************************************\n");
}


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
			fprintf(stderr, "********************************************************************************\n\n");
			fprintf(stderr, "Error: Could not open run-time I/O config file %s\n\n", fname);
			fprintf(stderr, "********************************************************************************\n");
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
 *  Function: attribute_check
 *
 *  Checks that a stream has the required attributes, and that attributes 
 *  are consistent.
 *
 *********************************************************************************/
int attribute_check(ezxml_t stream)
{
	const char *s_name, *s_type, *s_filename, *s_records, *s_input, *s_output;
	char msgbuf[MSGSIZE];
	int i, len, nextchar;


	s_name = ezxml_attr(stream, "name");
	s_type = ezxml_attr(stream, "type");
	s_filename = ezxml_attr(stream, "filename_template");
	s_records = ezxml_attr(stream, "records_per_file");
	s_input = ezxml_attr(stream, "input_interval");
	s_output = ezxml_attr(stream, "output_interval");


	/*
	 *  Check for required attributes
	 */
	if (s_name == NULL) {
		fmt_err("stream must have the \"name\" attribute.");
		return 1;
	}
	else if (s_type == NULL) {
		snprintf(msgbuf, MSGSIZE, "stream %s must have the \"type\" attribute.", s_name);
		fmt_err(msgbuf);
		return 1;
	}
	else if (s_filename == NULL) {
		snprintf(msgbuf, MSGSIZE, "stream %s must have the \"filename_template\" attribute.", s_name);
		fmt_err(msgbuf);
		return 1;
	}
	else if (s_records == NULL) {
		snprintf(msgbuf, MSGSIZE, "stream %s must have the \"records_per_file\" attribute.", s_name);
		fmt_err(msgbuf);
		return 1;
	}


	/*
	 *  Check that input streams have an input interval, output streams have an output interval
	 */
	if (strstr(s_type, "input") != NULL && s_input == NULL) {
		snprintf(msgbuf, MSGSIZE, "stream %s is an input stream and must have the \"input_interval\" attribute.", s_name);
		fmt_err(msgbuf);
		return 1;
	}
	else if (strstr(s_type, "output") != NULL && s_output == NULL) {
		snprintf(msgbuf, MSGSIZE, "stream %s is an output stream and must have the \"output_interval\" attribute.", s_name);
		fmt_err(msgbuf);
		return 1;
	}


	/*
	 *  Check that the filename template contains no illegal characters or variables
	 */
	len = strlen(s_filename);
	nextchar = 0;
	for (i=(len-1); i>=0; nextchar=s_filename[i--]) {
		if (s_filename[i] == '$') {
			if (strchr("YMDdhmsG",nextchar) == NULL) {
				snprintf(msgbuf, MSGSIZE, "filename_template for stream %s contains unrecognized variable \"$%c\".", s_name, nextchar);
				fmt_err(msgbuf);
				return 1;
			}
		}	
	}

	return 0;
}


/*********************************************************************************
 *
 *  Function: uniqueness_check
 *
 *  Checks that two streams have unique name and filename_template attributes
 *
 *********************************************************************************/
int uniqueness_check(ezxml_t stream1, ezxml_t stream2)
{
	const char *name, *name2;
	const char *filename, *filename2;
	char msgbuf[MSGSIZE];

	if (stream1 != stream2) {
		name = ezxml_attr(stream1, "name");
		filename = ezxml_attr(stream1, "filename_template");
		name2 = ezxml_attr(stream2, "name");
		filename2 = ezxml_attr(stream2, "filename_template");

		if (strcmp(name, name2) == 0) {
			snprintf(msgbuf, MSGSIZE, "stream \"%s\" is define more than once.", name);
			fmt_err(msgbuf);
			return 1;
		}
		if (strcmp(filename, filename2) == 0) {
			snprintf(msgbuf, MSGSIZE, "streams \"%s\" and \"%s\" cannot share the filename_template \"%s\".", name, name2, filename);
			fmt_err(msgbuf);
			return 1;
		}
	}

	return 0;
}


/*********************************************************************************
 *
 *  Function: check_streams
 *
 *  Validates the specification of run-time streams.
 *
 *********************************************************************************/
int check_streams(ezxml_t streams)
{
	ezxml_t stream_xml;
	ezxml_t stream2_xml;
	ezxml_t test_xml;
	ezxml_t test2_xml;
	const char *name;
	const char *filename;
	char msgbuf[MSGSIZE];


	/* Check immutable streams */
	for (stream_xml = ezxml_child(streams, "immutable_stream"); stream_xml; stream_xml = ezxml_next(stream_xml)) {
		if (attribute_check(stream_xml) != 0) {
			return 1;
		}	

		/* Check that users are not attempting to add fields to an immutable stream */
		test_xml = ezxml_child(stream_xml, "var");
		test2_xml = ezxml_child(stream_xml, "file");
		if (test_xml != NULL || test2_xml != NULL) {
			name = ezxml_attr(stream_xml, "name");
			snprintf(msgbuf, MSGSIZE, "the set of variables in stream \"%s\" cannot be modified.", name);
			fmt_err(msgbuf);
			return 1;
		}
	}

	/* Check mutable streams */
	for (stream_xml = ezxml_child(streams, "stream"); stream_xml; stream_xml = ezxml_next(stream_xml)) {
		name = ezxml_attr(stream_xml, "name");

		if (attribute_check(stream_xml) != 0) {
			return 1;
		}	

		/* If fields are specified in a separate file, that file should exist */
		for (test_xml = ezxml_child(stream_xml, "file"); test_xml; test_xml = ezxml_next(test_xml)) {
			filename = ezxml_attr(test_xml, "name");
			if (access(filename, F_OK|R_OK) == -1) {
				snprintf(msgbuf, MSGSIZE, "definition of stream \"%s\" references file %s that cannot be opened for reading.", name, filename);
				fmt_err(msgbuf);
				return 1;
			}
		}
	}


	/* Check that the name and filename_template attributes of all streams are unique */
	for (stream_xml = ezxml_child(streams, "stream"); stream_xml; stream_xml = ezxml_next(stream_xml)) {
		for (stream2_xml = ezxml_child(streams, "stream"); stream2_xml; stream2_xml = ezxml_next(stream2_xml)) {
			if (uniqueness_check(stream_xml, stream2_xml)) return 1;
		}
		for (stream2_xml = ezxml_child(streams, "immutable_stream"); stream2_xml; stream2_xml = ezxml_next(stream2_xml)) {
			if (uniqueness_check(stream_xml, stream2_xml)) return 1;
		}
	}
	for (stream_xml = ezxml_child(streams, "immutable_stream"); stream_xml; stream_xml = ezxml_next(stream_xml)) {
		for (stream2_xml = ezxml_child(streams, "stream"); stream2_xml; stream2_xml = ezxml_next(stream2_xml)) {
			if (uniqueness_check(stream_xml, stream2_xml)) return 1;
		}
		for (stream2_xml = ezxml_child(streams, "immutable_stream"); stream2_xml; stream2_xml = ezxml_next(stream2_xml)) {
			if (uniqueness_check(stream_xml, stream2_xml)) return 1;
		}
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
void xml_stream_parser(char *fname, void *manager, int *mpi_comm, int *status)
{
        char *xml_buf;
	size_t bufsize;
	ezxml_t streams;
	ezxml_t foo_xml;
	const char *attr;
	int err;

	fprintf(stderr, "MGD DEV begin parsing run-time I/O from %s\n", fname);
	*status = 0;

	global_file = fname;
        if (par_read(fname, mpi_comm, &xml_buf, &bufsize) != 0) {
		*status = 1;
		return;
	}

	streams = ezxml_parse_str(xml_buf, bufsize);
	if (!streams) {
		fprintf(stderr, "********************************************************************************\n\n");
		fprintf(stderr, "Error: Problems encountered while parsing run-time I/O config file %s\n", fname);
		fprintf(stderr, "********************************************************************************\n\n");
		*status = 1;
		return;
	}	

	if (check_streams(streams) != 0) {
		*status = 1;
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
