// Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
// and the University Corporation for Atmospheric Research (UCAR).
//
// Unless noted otherwise source code is licensed under the BSD license.
// Additional copyright and license information can be found in the LICENSE file
// distributed with this code, or at http://mpas-dev.github.com/license.html
//

#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include "ezxml/ezxml.h"

#define MSGSIZE 128


/* 
 *  Interface routines for building streams at run-time; defined in mpas_stream_manager.F
 */
void stream_mgr_create_stream_c(void *, const char *, int *, const char *, char *, char *, int *, int *, int *);
void mpas_stream_mgr_add_field_c(void *, const char *, const char *, int *);
void stream_mgr_add_alarm_c(void *, const char *, const char *, const char *, const char *, int *);


/*
 *  Stack node type used for basic syntax checking of XML
 */
struct stacknode {
	int line;
	char name[MSGSIZE];
	struct stacknode *next;
};

struct stacknode *head = NULL;


/* 
 *  Global variables
 */
static char *global_file;


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
 *  Function: fmt_warn
 *
 *  Prints a warning message in a standard format.
 *
 *********************************************************************************/
void fmt_warn(const char *mesg)
{
	fprintf(stderr,"********************************************************************************\n");
	fprintf(stderr,"Warning: In file %s, %s\n", global_file, mesg);
	fprintf(stderr,"********************************************************************************\n");
}


/*********************************************************************************
 *
 *  Function: push_tag
 *
 *  Pushes a new node onto the stack.
 *
 *********************************************************************************/
void push_tag(struct stacknode *node)
{
	if (node != NULL) {
		node->next = head;
		head = node;
	}
}


/*********************************************************************************
 *
 *  Function: pop_tag
 *
 *  Pops a new node from the stack.
 *
 *********************************************************************************/
struct stacknode * pop_tag(void)
{
	struct stacknode *retval;

	retval = head;
	if (head != NULL) {
		head = head->next;
	}

	return retval;
}


/*********************************************************************************
 *
 *  Function: parse_xml_tag_name
 *
 *  Copies only the name of an XML tag from tag_buf into tag_name. For example,
 *  the name of the tag
 * 
 *      <stream name="floop" interval="06:00:00"/>
 *
 *  is the string "stream".
 *
 *********************************************************************************/
void parse_xml_tag_name(char *tag_buf, char *tag_name)
{
	size_t i;

	/* Assume that a name ends with a space or null character */
	i = 0;
	while (tag_buf[i] != ' ' && tag_buf[i] != '\0') {
		tag_name[i] = tag_buf[i];
		i++;
	}

	tag_name[i] = '\0';
}


/*********************************************************************************
 *
 *  Function: parse_xml_tag
 *
 *  Parses the next XML tag into a string, plus other bookkeeping. All characters 
 *  between the first '<' and immediately following '>' character from xml_buf
 *  are copied into tag. The length of the buffer xml_buf is at least buf_len, and
 *  the length of the buffer tag is also at least buf_len.
 *
 *  For providing useful error messages, this routine also counts line numbers, 
 *  incrementing the line number each time a newline character is encountered.
 *
 *  The output argument start_line provides the line number on which the returned
 *  tag began.
 *
 *  The output argument tag_len provides the number of characters in the tag that
 *  were copied into the tag buffer. If no complete XML tag is found in the input 
 *  buffer, the tag_len argument will be set to 0.
 *
 *********************************************************************************/
size_t parse_xml_tag(char *xml_buf, size_t buf_len, char *tag, size_t *tag_len, int *line, int *start_line)
{
	size_t i, j;

	/* Look for beginning of tag */
	i = 0;
	while (i < buf_len && xml_buf[i] != '<') {
		if (xml_buf[i] == '\n')
			(*line)++;
		i++;
	}

	/* Ran out of characters... */
	if (i == buf_len) {
		*tag_len = 0;
		return 0;
	}


	/* Move on to next character after opening '<' */
	*start_line = *line;
	i++;

	/* Copy tag into string */
	j = 0;
	while (i < buf_len && xml_buf[i] != '>') {
		if (xml_buf[i] == '\n')
			(*line)++;
		tag[j] = xml_buf[i];
		i++;
		j++;
	}

	/* Didn't find a closing '>' character */
	if (i == buf_len) {
		*tag_len = 0;
		return 0;
	}

	tag[j] = '\0';
	i++;

	*tag_len = j;

	return i;
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
		snprintf(msgbuf, MSGSIZE, "stream \"%s\" must have the \"type\" attribute.", s_name);
		fmt_err(msgbuf);
		return 1;
	}
	else if (s_filename == NULL) {
		snprintf(msgbuf, MSGSIZE, "stream \"%s\" must have the \"filename_template\" attribute.", s_name);
		fmt_err(msgbuf);
		return 1;
	}
	else if (s_records == NULL) {
		snprintf(msgbuf, MSGSIZE, "stream \"%s\" must have the \"records_per_file\" attribute.", s_name);
		fmt_err(msgbuf);
		return 1;
	}


	/*
	 *  Check that input streams have an input interval, output streams have an output interval
	 */
	if (strstr(s_type, "input") != NULL && s_input == NULL) {
		snprintf(msgbuf, MSGSIZE, "stream \"%s\" is an input stream and must have the \"input_interval\" attribute.", s_name);
		fmt_err(msgbuf);
		return 1;
	}
	if (strstr(s_type, "output") != NULL && s_output == NULL) {
		snprintf(msgbuf, MSGSIZE, "stream \"%s\" is an output stream and must have the \"output_interval\" attribute.", s_name);
		fmt_err(msgbuf);
		return 1;
	}
	if (strstr(s_type, "input") != NULL && strstr(s_type, "output") == NULL && s_output != NULL) {
		snprintf(msgbuf, MSGSIZE, "input-only stream \"%s\" has the \"output_interval\" attribute.", s_name);
		fmt_warn(msgbuf);
	}
	if (strstr(s_type, "output") != NULL && strstr(s_type, "input") == NULL && s_input != NULL) {
		snprintf(msgbuf, MSGSIZE, "output-only stream \"%s\" has the \"input_interval\" attribute.", s_name);
		fmt_warn(msgbuf);
	}


	/*
	 *  Check that the filename template contains no illegal characters or variables
	 *  NB: If new variable characters are added here, they should also be accommodated in
	 *      the mpas_expand_string() subroutine in the mpas_timekeeping module.
	 */
	len = strlen(s_filename);
	nextchar = 0;
	for (i=(len-1); i>=0; nextchar=s_filename[i--]) {
		if (s_filename[i] == '$') {
			if (strchr("YMDdhmsG",nextchar) == NULL) {
				snprintf(msgbuf, MSGSIZE, "filename_template for stream \"%s\" contains unrecognized variable \"$%c\".", s_name, nextchar);
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
/* TODO: should this also be done only on the master task? */
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
 *  Function: xml_syntax_check
 *
 *  Performs a few basic syntax checks on a buffer containing XML:
 *  1) Are the angle brackets balanced?
 *  2) Are the quotes balanced?
 *  3) Are all XML tags closed and nested properly?
 *
 *  There are clearly many syntax errors that this code will not catch, e.g., attribute
 *  values that contain no quotes at all; similarly, there are syntactically correct
 *  situations that this code will flag as bad, e.g., quoted strings that contain the
 *  '=' character. If we really wanted to be thorough, we should employ a proper parser
 *  with a well-specified grammar.
 *
 *********************************************************************************/
int xml_syntax_check(char *xml_buf, size_t bufsize)
{
	size_t i;
	size_t len;
	int nleft, nright, line, start_line;
	char msgbuf[MSGSIZE];
	char *tag_buf;
	struct stacknode *node;
	struct stacknode tmp_node;


	/* 
	 *  Check that we have balanced angle brackets 
	 */
	nleft = 0;
	nright = 0;
	line = 1;
	for (i=0; i<bufsize; i++) {
		if (xml_buf[i] == '<') {
			nleft++;
			if (nleft - nright > 1) {
				snprintf(msgbuf, MSGSIZE, "line %i, unexpected \'<\' character. Is the previous XML tag missing a \'>\'?", line);
				fmt_err(msgbuf);
				return 1;
			}
		}
		else if (xml_buf[i] == '>') {
			nright++;
			if (nleft != nright) {
				snprintf(msgbuf, MSGSIZE, "line %i, unexpected \'>\' character. Is the XML tag missing a \'<\'?", line);
				fmt_err(msgbuf);
				return 1;
			}
		}
		else if (xml_buf[i] == '\n') {
			line++;
		}
	}
	if (nleft != nright) {                                  /* Probably only triggered if no final '>' character? */
		fmt_err("unbalanced angle brackets in XML. Is the file missing a final \'>\'?");
		return 1;
	}


	/*
	 *  Check that we have balanced quotes
	 *  NB: This simple logic WILL NOT WORK if quoted strings are allowed to contain the '=' character!
	 */
	nleft = 0;
	line = 1;
	for (i=0; i<bufsize; i++) {
		if (xml_buf[i] == '"') {
			nleft = (nleft + 1) % 2;
		}

		/* 
		 *  When we reach the end of a line or the beginning of a new attribute definition,
		 *     the quotes should be balanced...
		 */
		if (xml_buf[i] == '=' || xml_buf[i] == '\n') {
			if (nleft != 0) {
				snprintf(msgbuf, MSGSIZE, "line %i, unterminated string. Is a closing quote not present?", line);
				fmt_err(msgbuf);
				return 1;
			}
			if (xml_buf[i] == '\n') {
				line++;
			}
		}
	}
	if (nleft != 0) {
		fmt_err("unbalanced quotes in XML.");
		return 1;
	}


	/* 
	 *  Check that each tag is closed
	 */
	i = 0;
	tag_buf = (char *)malloc(bufsize);

	line = 1;
	do {
		i += parse_xml_tag(&xml_buf[i], bufsize, tag_buf, &len, &line, &start_line);

		if (len > 0) {

			/* Probably a comment tag -- though this is not a perfect check... */
			if (tag_buf[0] == '!' && tag_buf[len-1] == '-') {
				/* Would it be better for the tag parser to just skip over comments? */
			}
			/* An opening tag. Push it onto the stack... */
			else if (tag_buf[0] != '/' && tag_buf[len-1] != '/') {
				node = (struct stacknode *)malloc(sizeof(struct stacknode));
				parse_xml_tag_name(tag_buf, node->name);
				node->line = start_line;
				push_tag(node);
			}
			/* A closing tag. Pop the stack... */
			else if (tag_buf[0] == '/' && tag_buf[len-1] != '/') {
				node = pop_tag();
				parse_xml_tag_name(&tag_buf[1], tmp_node.name);    /* NB: &tag_buf[1] to skip over '/' character */
				if (strncmp(tmp_node.name, node->name, (size_t)MSGSIZE) != 0) {
					fprintf(stderr, "Found unexpected closing tag \"%s\" at line %i.\n", tmp_node.name, start_line);
					snprintf(msgbuf, MSGSIZE, "line %i, unclosed or badly nested XML tag \"%s\".", node->line, node->name);
					fmt_err(msgbuf);

					while ((node = pop_tag()) != NULL) 
						free(node);
					return 1;	
				}
				free(node);
			}
			/* A singleton tag. Life is simple... */
			else if (tag_buf[0] != '/' && tag_buf[len-1] == '/') {
				parse_xml_tag_name(tag_buf, tmp_node.name);

			}
			/* Probable syntax error? */
			else {
				
			}

		}

	} while (i <= bufsize && len > 0);

	/* Pop the rest of the stack for any unclosed tags */
	node = pop_tag();
	if (node != NULL) {
		snprintf(msgbuf, MSGSIZE, "line %i, unclosed or badly nested XML tag \"%s\".", node->line, node->name);
		fmt_err(msgbuf);
		
		while ((node = pop_tag()) != NULL) 
			free(node);
		return 1;	
	}

	free(tag_buf);

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
	ezxml_t stream_xml;
	ezxml_t varfile_xml;
	ezxml_t var_xml;
	const char *streamID, *filename_template, *direction, *records, *varfile, *fieldname_const, *reference_time, *record_interval;
	const char *interval_in, *interval_out;
	char ref_time_local[256];
	char rec_intv_local[256];
	char fieldname[256];
	FILE *fd;
	char msgbuf[MSGSIZE];
	int itype;
	int irecs;
	int immutable;
	int err;


	fprintf(stderr, "\nParsing run-time I/O configuration from %s ...\n", fname);
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

	err = 0;

	/* First, handle changes to immutable stream filename templates, intervals, etc. */
	immutable = 1;
	for (stream_xml = ezxml_child(streams, "immutable_stream"); stream_xml; stream_xml = ezxml_next(stream_xml)) {
		streamID = ezxml_attr(stream_xml, "name");
		direction = ezxml_attr(stream_xml, "type");
		filename_template = ezxml_attr(stream_xml, "filename_template");
		records = ezxml_attr(stream_xml, "records_per_file");
		interval_in = ezxml_attr(stream_xml, "input_interval");
		interval_out = ezxml_attr(stream_xml, "output_interval");
		reference_time = ezxml_attr(stream_xml, "reference_time");
		record_interval = ezxml_attr(stream_xml, "record_interval");

		fprintf(stderr, "   found immutable stream %s\n", streamID);

		irecs = atoi(records);

		/* NB: These type constants must match those in the mpas_stream_manager module! */
		if (strstr(direction, "input") != NULL && strstr(direction, "output") != NULL)
			itype = 3;
		else if (strstr(direction, "input") != NULL)
			itype = 1;
		else if (strstr(direction, "output") != NULL) 
			itype = 2;
		else 
			itype = 4;

		if (reference_time != NULL) {
			snprintf(ref_time_local, 256, "%s", reference_time);
		}
		else {
			snprintf(ref_time_local, 256, "none");
		}

		if (record_interval != NULL) {
			snprintf(rec_intv_local, 256, "%s", record_interval);
		}
		else {
			snprintf(rec_intv_local, 256, "none");
		}

		stream_mgr_create_stream_c(manager, streamID, &itype, filename_template, ref_time_local, rec_intv_local, 
                                           &irecs, &immutable, &err);
		if (err != 0) {
			*status = 1;
			return;
		}

		/* Possibly add an input alarm for this stream */
		if (itype == 3 || itype == 1) {
			stream_mgr_add_alarm_c(manager, streamID, "input", "start", interval_in, &err);
			if (err != 0) {
				*status = 1;
				return;
			}
		}

		/* Possibly add an output alarm for this stream */
		if (itype == 3 || itype == 2) {
			stream_mgr_add_alarm_c(manager, streamID, "output", "start", interval_out, &err);
			if (err != 0) {
				*status = 1;
				return;
			}
		}
	}

	/* Next, handle modifications to mutable streams as well as new stream definitions */
	immutable = 0;
	for (stream_xml = ezxml_child(streams, "stream"); stream_xml; stream_xml = ezxml_next(stream_xml)) {
		streamID = ezxml_attr(stream_xml, "name");
		direction = ezxml_attr(stream_xml, "type");
		filename_template = ezxml_attr(stream_xml, "filename_template");
		records = ezxml_attr(stream_xml, "records_per_file");
		interval_in = ezxml_attr(stream_xml, "input_interval");
		interval_out = ezxml_attr(stream_xml, "output_interval");
		reference_time = ezxml_attr(stream_xml, "reference_time");
		record_interval = ezxml_attr(stream_xml, "record_interval");

		fprintf(stderr, "   found mutable stream %s\n", streamID);

		irecs = atoi(records);

		/* NB: These type constants must match those in the mpas_stream_manager module! */
		if (strstr(direction, "input") != NULL && strstr(direction, "output") != NULL)
			itype = 3;
		else if (strstr(direction, "input") != NULL)
			itype = 1;
		else if (strstr(direction, "output") != NULL) 
			itype = 2;
		else 
			itype = 4;

		if (reference_time != NULL) {
			snprintf(ref_time_local, 256, "%s", reference_time);
		}
		else {
			snprintf(ref_time_local, 256, "none");
		}

		if (record_interval != NULL) {
			snprintf(rec_intv_local, 256, "%s", record_interval);
		}
		else {
			snprintf(rec_intv_local, 256, "none");
		}


		stream_mgr_create_stream_c(manager, streamID, &itype, filename_template, ref_time_local, rec_intv_local, 
                                           &irecs, &immutable, &err);
		if (err != 0) {
			*status = 1;
			return;
		}

		/* Possibly add an input alarm for this stream */
		if (itype == 3 || itype == 1) {
			stream_mgr_add_alarm_c(manager, streamID, "input", "start", interval_in, &err);
			if (err != 0) {
				*status = 1;
				return;
			}
		}

		/* Possibly add an output alarm for this stream */
		if (itype == 3 || itype == 2) {
			stream_mgr_add_alarm_c(manager, streamID, "output", "start", interval_out, &err);
			if (err != 0) {
				*status = 1;
				return;
			}
		}

		for (varfile_xml = ezxml_child(stream_xml, "file"); varfile_xml; varfile_xml = ezxml_next(varfile_xml)) {
			varfile = ezxml_attr(varfile_xml, "name");
/* TODO: We should probably only have one task open and read the file... */
			fd = fopen(varfile, "r");
			if (fd != NULL) {
				while (fscanf(fd, "%s", fieldname) != EOF) {
					stream_mgr_add_field_c(manager, streamID, (const char *)fieldname, &err);
					if (err != 0) {
						*status = 1;
						return;
					}
				}
				fclose(fd);
			}
			else {
				snprintf(msgbuf, MSGSIZE, "definition of stream \"%s\" references file %s that cannot be opened for reading.", streamID, varfile);
				fmt_err(msgbuf);
				*status = 1;
				return;
			}
		}
		for (var_xml = ezxml_child(stream_xml, "var"); var_xml; var_xml = ezxml_next(var_xml)) {
			fieldname_const = ezxml_attr(var_xml, "name");
			stream_mgr_add_field_c(manager, streamID, fieldname_const, &err);
			if (err != 0) {
				*status = 1;
				return;
			}
		}
	}

	free(xml_buf);

	fprintf(stderr, "... done parsing run-time I/O\n\n");
}


/*********************************************************************************
 *
 *  Function: xml_stream_get_filename
 *
 *  Parses an XML file and searches for the stream whose name matches the 'streamname'
 *  argument; then, returns the associated filename template for that stream in 
 *  the 'filename' argument.
 *
 *  The fname argument provides the name of the XML file that contains the stream
 *  definitions, and mpi_comm is the Fortran MPI communicator used by MPAS.
 *
 *********************************************************************************/
void xml_stream_get_filename(char *fname, char *streamname, int *mpi_comm, char *filename, int *status)
{
	char *xml_buf;
	size_t bufsize;
	ezxml_t streams;
	ezxml_t stream_xml;
	const char *streamID, *filename_template;
	int found;


	*status = 0;

	global_file = fname;
	if (par_read(fname, mpi_comm, &xml_buf, &bufsize) != 0) {
		*status = 1;
		return;
	}

	if (xml_syntax_check(xml_buf, bufsize) != 0) {
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

	found = 0;
	for (stream_xml = ezxml_child(streams, "immutable_stream"); stream_xml; stream_xml = ezxml_next(stream_xml)) {
		streamID = ezxml_attr(stream_xml, "name");
		filename_template = ezxml_attr(stream_xml, "filename_template");

		if (strcmp(streamID, streamname) == 0) {
			found = 1;
			fprintf(stderr, "Found grid stream with template %s\n", filename_template);
			sprintf(filename, "%s", filename_template);
			break;
		}
	}
	if (found == 0) {
		*status = 1;
		return;
	}
}
