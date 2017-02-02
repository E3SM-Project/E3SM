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
#include <errno.h>
#include "ezxml.h"

#ifdef _MPI
#include "mpi.h"
#endif

#define MSGSIZE 256


/*
 *  Interface routines for building streams at run-time; defined in mpas_stream_manager.F
 */
void stream_mgr_create_stream_c(void *, const char *, int *, const char *, const char *, const char *, const char *, int *, int *, int *, int *, int *);
void mpas_stream_mgr_add_field_c(void *, const char *, const char *, const char *, int *);
void mpas_stream_mgr_add_immutable_stream_fields_c(void *, const char *, const char *, const char *, int *);
void mpas_stream_mgr_add_pool_c(void *, const char *, const char *, const char *, int *);
void stream_mgr_add_alarm_c(void *, const char *, const char *, const char *, const char *, int *);
void stream_mgr_add_pkg_c(void *, const char *, const char *, int *);


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
	fprintf(stderr,"* Error: In file %s, %s\n", global_file, mesg);
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
	fprintf(stderr,"* Warning: In file %s, %s\n", global_file, mesg);
	fprintf(stderr,"********************************************************************************\n");
}

/*********************************************************************************
 *
 *  Function: fmt_info
 *
 *  Prints an informational message in a standard format.
 *
 *********************************************************************************/
void fmt_info(const char *mesg)
{
	fprintf(stderr,"\n");
	fprintf(stderr," Information: In file %s, %s\n", global_file, mesg);
	fprintf(stderr,"\n");
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
 *  The return value is the index in xml_buf representing the end of the tag,
 *  realtive to the starting position.
 *
 *********************************************************************************/
size_t parse_xml_tag(char *xml_buf, size_t buf_len, char *tag, size_t *tag_len, int *line, int *start_line)
{
	size_t i, j;
	int found_end, block_comment;


	i = 0;
	do {
		/* Look for beginning of tag */
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

		block_comment = 0;
		/* Skip comment tags */
		if ( xml_buf[i] == '!' && xml_buf[i+1] == '-' && xml_buf[i+2] == '-' ) {
			block_comment = 1;

			/* find end of the comment... */
			i = i+2;
			found_end = 0;
			while (i < buf_len && ! found_end) {
				if ( xml_buf[i] == '-' && xml_buf[i+1] == '-' && xml_buf[i+2] == '>' ) {
					found_end = 1;
					i = i+2;
				} else if ( xml_buf[i] == '\n' ) {
					(*line)++;
				}

				i++;
			}


			/* Ran out of characters... */
			if (i == buf_len) {
				*tag_len = 0;
				return 0;
			}
		}
	} while (block_comment);

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
	MPI_Comm comm;

	comm = MPI_Comm_f2c((MPI_Fint)(*mpi_comm));
	err = MPI_Comm_rank(comm, &rank);
#else
	rank = 0;
#endif

	if (rank == 0) {
		iofd = open(fname, O_RDONLY);
		if (iofd <= 0) {
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
		err = read(iofd, (void *)(*xml_buf), *bufsize);

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
	const char *s_name, *s_type, *s_filename, *s_filename_intv, *s_input, *s_output, *s_ref_time;
	char msgbuf[MSGSIZE];
	int i, len, nextchar;

	s_name = ezxml_attr(stream, "name");
	s_type = ezxml_attr(stream, "type");
	s_filename = ezxml_attr(stream, "filename_template");
	s_filename_intv = ezxml_attr(stream, "filename_interval");
	s_input = ezxml_attr(stream, "input_interval");
	s_output = ezxml_attr(stream, "output_interval");
	s_ref_time = ezxml_attr(stream, "reference_time");


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
	 *  Check that filename_interval is given an acceptable value.
	 */
	if ( s_filename_intv != NULL ) {
		if ( strstr(s_filename_intv, "input_interval") != NULL && s_input == NULL) {
			snprintf(msgbuf, MSGSIZE, "stream \"%s\" has a value of \"input_interval\" for the \"filename_interval\" attribute, without defining the \"input_interval\" attribute.", s_name);
			fmt_err(msgbuf);
			return 1;
		}
		if ( strstr(s_filename_intv, "output_interval") != NULL && s_output == NULL) {
			snprintf(msgbuf, MSGSIZE, "stream \"%s\" has a value of \"output_interval\" for the \"filename_interval\" attribute, without defining the \"output_interval\" attribute.", s_name);
			fmt_err(msgbuf);
			return 1;
		}
		if ( strstr(s_filename_intv, "input_interval") != NULL && strstr(s_input, "initial_only") != NULL) {
			snprintf(msgbuf, MSGSIZE, "stream \"%s\" cannot have a value of \"input_interval\" for the \"filename_interval\" attribute, when \"input_interval\" is set to \"initial_only\".", s_name);
			fmt_err(msgbuf);
			return 1;
		}
		if ( strstr(s_filename_intv, "output_interval") != NULL && strstr(s_output, "initial_only") != NULL) {
			snprintf(msgbuf, MSGSIZE, "stream \"%s\" cannot have a value of \"output_interval\" for the \"filename_interval\" attribute, when \"output_interval\" is set to \"initial_only\".", s_name);
			fmt_err(msgbuf);
			return 1;
		}
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
			if (strchr("YMDdhmsGSB",nextchar) == NULL) {
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
	const char *type, *type2;
	char msgbuf[MSGSIZE];

	if (stream1 != stream2) {
		name = ezxml_attr(stream1, "name");
		filename = ezxml_attr(stream1, "filename_template");
		type = ezxml_attr(stream1, "type");
		name2 = ezxml_attr(stream2, "name");
		filename2 = ezxml_attr(stream2, "filename_template");
		type2 = ezxml_attr(stream2, "type");

		if (strcmp(name, name2) == 0) {
			snprintf(msgbuf, MSGSIZE, "stream \"%s\" is define more than once.", name);
			fmt_err(msgbuf);
			return 1;
		}
		if (strstr(type, "output") != NULL || strstr(type2, "output") != NULL){
			if (strcmp(filename, filename2) == 0) {
				snprintf(msgbuf, MSGSIZE, "Output streams \"%s\" and \"%s\" cannot share the filename_template \"%s\".", name, name2, filename);
				fmt_err(msgbuf);
				return 1;
			}
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
	int nleftcom, nrightcom;
	char msgbuf[MSGSIZE];
	char *tag_buf;
	struct stacknode *node;
	struct stacknode tmp_node;


	/*
	 *  Check that we have balanced angle brackets
	 */
	nleft = 0;
	nright = 0;
	nleftcom = 0;
	nrightcom = 0;
	line = 1;

	if ( xml_buf[0] == '>' ) {
		snprintf(msgbuf, MSGSIZE, "line %i, unexpected starting \'>\' character. A file cannot start with a  \'>\' character.", line);
		fmt_err(msgbuf);
		return 1;
	}

	for (i=0; i<bufsize; i++) {
		if (xml_buf[i] == '<') {
			if (i+1 < bufsize && xml_buf[i+1] == '!'){
				nleftcom++;
				if (nleftcom - nrightcom > 1) {
					snprintf(msgbuf, MSGSIZE, "line %i, unexpected XML comment open. Is the previous XML comment missing a \'-->\'?", line);
					fmt_err(msgbuf);
					return 1;
				} else if (nleft != nright) {
					snprintf(msgbuf, MSGSIZE, "line %i, unexpected XML comment open. Is the previous XML tag missing a \'>\'?\n   NOTE: Comments are not allowed within an open XML tag.", line);
					fmt_err(msgbuf);
					return 1;
				}
			} else {
				nleft++;
				if (nleft - nright > 1){
					snprintf(msgbuf, MSGSIZE, "line %i, unexpected \'<\' character. Is the previous XML tag missing a \'>\'?", line);
					fmt_err(msgbuf);
					return 1;
				}
			}
		}
		else if (xml_buf[i] == '>') {
			if (i > 0 && xml_buf[i-1] == '-'){
				nrightcom++;
				if (nleftcom != nrightcom) {
					snprintf(msgbuf, MSGSIZE, "line %i, unexpected XML comment close. Is the XML comment missing a \'<!--\'?", line);
					fmt_err(msgbuf);
					return 1;
				} else if (nleft != nright) {
					snprintf(msgbuf, MSGSIZE, "line %i, unexpected XML comment close. Is the previous XML tag missing a \'>\'?\n   NOTE: Comments are not allowed within an open XML tag.", line);
					fmt_err(msgbuf);
					return 1;
				}
			} else {
				nright++;
				if (nleft != nright) {
					snprintf(msgbuf, MSGSIZE, "line %i, unexpected \'>\' character. Is the XML tag missing a \'<\'?", line);
					fmt_err(msgbuf);
					return 1;
				}
			}
		}
		else if (xml_buf[i] == '\n') {
			line++;
		}
	}
	if (nleft != nright) {				  /* Probably only triggered if no final '>' character? */
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
		i += parse_xml_tag(&xml_buf[i], (bufsize - i), tag_buf, &len, &line, &start_line);

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
 *  Function: build_stream_path
 *
 *  Takes as input a string defining the filename template for a stream, and, for
 *  each directory in the template, ensures that the directory exists. If a directory
 *  in the template already exists but is not writable, a non-zero error is returned.
 *
 *********************************************************************************/
int build_stream_path(const char *stream, const char *template, int *mpi_comm)
{
	char *filename_path;
	char *directory;
	int create_dir;
	int i, len;
	char msgbuf[MSGSIZE];
	int err, retval;
	int writable_parent;
	int rank;


#ifdef _MPI
	MPI_Comm comm;

	comm = MPI_Comm_f2c((MPI_Fint)(*mpi_comm));
	err = MPI_Comm_rank(comm, &rank);
#else
	rank = 0;
#endif

	if (rank == 0) {
		/*
		 * Check that paths in a filename template exist for output streams.
		 * Create them if they don't.
		 *
		 * Parse immutable steams first
		 */
		create_dir = 0;
		filename_path = strdup(template);

		len = strlen(filename_path);
		directory = (char *)malloc(sizeof(char) * (size_t)len);

		for (i=(len-1); i>=0; i--) {
			if (filename_path[i] == '/') {
				filename_path[i] = '\0';
				create_dir = 1;
				break;
			}
		}

		if (create_dir) {
			writable_parent = 1;
			for(i=0; i < len; i++) {
				if (filename_path[i] == '/') {
					directory[i] = '\0';
					err = mkdir(directory, S_IRWXU | S_IRWXG | S_IRWXO);
					if ( err != 0 ) {
						if ( errno == EEXIST ) {
							/* directory exists, need to check permissions */
							writable_parent = 1;
							if (access(directory, W_OK) != 0) {
								writable_parent = 0;
							}
						} else if ( !writable_parent ) {
							snprintf(msgbuf, MSGSIZE, "cannot create directory %s needed by stream \"%s\": parent directory is not writable.", directory, stream);
							fmt_err(msgbuf);
							free(filename_path);
							free(directory);
							writable_parent = 0;

							retval = 1;
#ifdef _MPI
							err = MPI_Bcast(&retval, 1, MPI_INT, 0, comm);
#endif
							return retval;						
						}
					}
				}
				directory[i] = filename_path[i];
			}
			err = mkdir(filename_path, S_IRWXU | S_IRWXG | S_IRWXO);
			if ( ! err ) {
				fprintf(stderr, "        *** created directory %s for stream \"%s\"\n", filename_path, stream);
			} else if ( errno == EEXIST ) {
				/* directory exists, need to check permissions */
				if (access(filename_path, W_OK) != 0) {
					snprintf(msgbuf, MSGSIZE, "definition of stream \"%s\" references directory %s without write permission.", stream, filename_path);
					fmt_err(msgbuf);
					free(filename_path);
					free(directory);
	
					retval = 1;
#ifdef _MPI
					err = MPI_Bcast(&retval, 1, MPI_INT, 0, comm);
#endif
					return retval;
				}
			}
			else if ( !writable_parent ) {
					snprintf(msgbuf, MSGSIZE, "cannot create directory %s needed by stream \"%s\": parent directory is not writable.", directory, stream);
					fmt_err(msgbuf);
					free(filename_path);
					free(directory);
					writable_parent = 0;

					retval = 1;
#ifdef _MPI
					err = MPI_Bcast(&retval, 1, MPI_INT, 0, comm);
#endif
					return retval;						
			}
		}

		free(filename_path);
		free(directory);

		retval = 0;
#ifdef _MPI
		err = MPI_Bcast(&retval, 1, MPI_INT, 0, comm);
#endif
	}
#ifdef _MPI
	else {
		err = MPI_Bcast(&retval, 1, MPI_INT, 0, comm);
		if (retval != 0) {
			fprintf(stderr, "********************************************************************************\n");
			fprintf(stderr, "* Please check the standard error log from task 0 for error messages.\n");
			fprintf(stderr, "********************************************************************************\n");
		}
	}
#endif

	return retval;
}


/*********************************************************************************
 *
 *  Function: extract_stream_interval
 *
 *  Given an interval specification for a stream (interval) that references 
 *  an interval in another stream (e.g., "stream:history:output_interval"), and 
 *  an interval type (interval_type, either "input_interval" or "output_interval"), 
 *  extracts the value of the interval from the other stream and returns it in 
 *  the output argument interval2.
 *
 *  If the interval specification in the interval argument does not reference 
 *  another stream, the contents of interval2 are unchanged upon return from 
 *  this function.
 *
 *  In case the input interval references an interval in another stream and this 
 *  interval cannot found, this function returns a value of 1; otherwise, this 
 *  function returns 0.
 *
 *********************************************************************************/
int extract_stream_interval(const char *interval, const char *interval_type, const char **interval2, const char *streamID, ezxml_t streams)
{
	int i;
	int stream_found, copy_start, copy_from, copy_to;
	char match_stream_name[256];
	char interval_name[256];
	const char *streamID2;
	ezxml_t stream2_xml;
	ezxml_t streammatch_xml;


	if ( strncmp(interval, "stream:", 7) == 0 ) {

		/* Extract the name of the stream, and the name of the interval to use for interval */
		snprintf(match_stream_name, 256, "%s", (interval)+7);
		copy_start = -1;
		copy_from = -1;
		copy_to = 0;
		for ( i = 0; i < strlen(match_stream_name); i++ ) {
			if ( match_stream_name[i] == ':' ) {
				copy_start = i;
				copy_from = copy_start+1;
			}

			if ( copy_from == i ) {
				interval_name[copy_to] = match_stream_name[copy_from];
				copy_from++;
				copy_to++;
			}
		}
		match_stream_name[copy_start] = '\0';
		interval_name[copy_to] = '\0';

		if ( strcmp(match_stream_name, streamID) == 0 && strcmp(interval_name, interval_type) == 0 ) {
			fprintf(stderr, "ERROR: Attribute '%s' of stream '%s' references itself.\n", interval_type, streamID);
			return 1;
		}

		if ( strcmp(interval_name, "input_interval") != 0 && strcmp(interval_name, "output_interval") != 0 ) {
			fprintf(stderr, "ERROR: Attribute '%s' of stream '%s' references an invalid attribute: '%s'.\n", interval_type, streamID, interval_name);
			fprintf(stderr, "       Valid attributes are 'input_interval' and 'output_interval'.\n");
			return 1;
		}

		stream_found = 0;
		for ( stream2_xml = ezxml_child(streams, "immutable_stream"); stream2_xml && !stream_found; stream2_xml = stream2_xml->next ) {
			streamID2 = ezxml_attr(stream2_xml, "name");

			if ( strcmp(streamID2, match_stream_name) == 0 ){
				stream_found = 1;
				streammatch_xml = stream2_xml;
			}
		}

		for ( stream2_xml = ezxml_child(streams, "stream"); stream2_xml && !stream_found; stream2_xml = stream2_xml->next ) {
			streamID2 = ezxml_attr(stream2_xml, "name");

			if ( strcmp(streamID2, match_stream_name) == 0 ) {
				stream_found = 1;
				streammatch_xml = stream2_xml;
			}
		}

		if ( stream_found == 1 ) {
			*interval2 = ezxml_attr(streammatch_xml, interval_name);
		} 
		else {
			fprintf(stderr, "ERROR: The '%s' attribute of stream '%s' refers to an undefined stream named '%s'.\n", interval_type, streamID, match_stream_name);
			return 1;
		}


		if ( *interval2 == NULL ) {
			fprintf(stderr, "ERROR: The '%s' attribute of stream '%s' refers to an undefined attribute named '%s' of stream '%s'.\n", interval_type, streamID, interval_name, match_stream_name);
			return 1;
		} 
		else if ( strcmp(*interval2, "input_interval") == 0 || strcmp(*interval2, "output_interval") == 0 || strncmp(*interval2, "stream:", 7) == 0 ) {
			fprintf(stderr, "ERROR: The '%s' attribute of stream '%s' contains an unexpandable value: '%s'.\n", interval_type, streamID, *interval2);
			return 1;
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
	ezxml_t stream_xml;
	ezxml_t varfile_xml;
	ezxml_t var_xml;
	ezxml_t vararray_xml;
	ezxml_t varstruct_xml;
	ezxml_t substream_xml;
	ezxml_t stream2_xml;
	ezxml_t streammatch_xml;
	const char *compstreamname_const, *structname_const;
	const char *streamID, *filename_template, *filename_interval, *direction, *varfile, *fieldname_const, *reference_time, *record_interval, *streamname_const, *precision;
	const char *interval_in, *interval_out, *packagelist;
	const char *clobber;
	const char *iotype;
	const char *streamID2, *interval_in2, *interval_out2;
	char interval_name[256];
	char match_stream_name[256];
	char *packages, *package;
	char filename_interval_string[256];
	char ref_time_local[256];
	char rec_intv_local[256];
	char fieldname[256];
	char packages_local[256];
	char interval_type[32];
	FILE *fd;
	char msgbuf[MSGSIZE];
	int itype;
	int iclobber;
	int i_iotype;
	int iprec;
	int immutable;
	int stream_found, copy_start, copy_from, copy_to;
	int i;
	int err;


	packages_local[0] = '\0';

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
		filename_interval = ezxml_attr(stream_xml, "filename_interval");
		interval_in = ezxml_attr(stream_xml, "input_interval");
		interval_in2 = ezxml_attr(stream_xml, "input_interval");
		interval_out = ezxml_attr(stream_xml, "output_interval");
		interval_out2 = ezxml_attr(stream_xml, "output_interval");
		reference_time = ezxml_attr(stream_xml, "reference_time");
		record_interval = ezxml_attr(stream_xml, "record_interval");
		precision = ezxml_attr(stream_xml, "precision");
		packagelist = ezxml_attr(stream_xml, "packages");
		clobber = ezxml_attr(stream_xml, "clobber_mode");
		iotype = ezxml_attr(stream_xml, "io_type");

		/* Extract the input interval, if it refer to other streams */
		if ( interval_in ) {
			sprintf(interval_type, "input_interval");
			*status = extract_stream_interval(interval_in, interval_type, &interval_in2, streamID, streams);
			if ( *status != 0 ) {
				return;
			}
		}

		/* Extract the output interval, if it refer to other streams */
		if ( interval_out ) {
			sprintf(interval_type, "output_interval");
			*status = extract_stream_interval(interval_out, interval_type, &interval_out2, streamID, streams);
			if ( *status != 0 ) {
				return;
			}
		}

		/* Setup filename_interval correctly.
		 *
		 * If filename_interval is not explicitly set...
		 *  - Default to input interval if it is set and an interval (for input, and input;output streams)
		 *  - Default to output interval if it is set and an interval (for output streams)
		 *  - Default to none for all other cases
		 *
		 * After this check is complete, if filename_interval still has a value of NULL, it will be replaced with 'none'
		 */
		if ( filename_interval == NULL){
			/* Check for an input;output stream. Handle first as this is the most complicated case. */
			if ( strstr(direction, "input") != NULL && strstr(direction, "output") != NULL ) {

				/* If input interval is an interval (i.e. not initial_only or none) set filename_interval to the interval. */
				if ( strstr(interval_in, "initial_only") == NULL && strstr(interval_in, "none") == NULL ){
					filename_interval = interval_in2;

				/* If output interval is an interval (i.e. not initial_only or none) set filename_interval to the interval. */
				} else if ( strstr(interval_out, "initial_only") == NULL && strstr(interval_out, "none") == NULL ){
					filename_interval = interval_out2;
				}
			/* Check for an input stream. */
			} else if ( strstr(direction, "input") != NULL ) {
				if ( strstr(interval_in2, "initial_only") == NULL && strstr(interval_in2, "none") == NULL ){
					filename_interval = interval_in2;
				}

			/* Check for an output stream. */
			} else if ( strstr(direction, "output") != NULL ) {
				if ( strstr(interval_out, "initial_only") == NULL && strstr(interval_out, "none") == NULL ){
					filename_interval = interval_out2;
				}
			}
		} else {
			/* Handle the case where filename_interval has a value of either:
			 * output_interval -- Overwrite filename_interval with the interval provided in output_interval
			 * input_interval -- Overwrite filename_interval with the interval provided in input_interval
			 *
			 * In either of these cases, if the intervals are set to be initial_only or none, nullify filename_interval
			 * to force it's value to be none as well.
			 */
			if ( strstr(filename_interval, "input_interval") != NULL ) {
				if ( strstr(interval_in, "initial_only") == NULL && strstr(interval_in, "none") == NULL ) {
					filename_interval = interval_in2;
				} else {
					filename_interval = NULL;
				}
			} else if ( strstr(filename_interval, "output_interval") != NULL ) {
				if ( strstr(interval_out, "initial_only") == NULL && strstr(interval_out, "none") == NULL ) {
					filename_interval = interval_out2;
				} else {
					filename_interval = NULL;
				}
			}
		}

		if ( filename_interval == NULL ) {
			sprintf(filename_interval_string, "none");
		} else {
			sprintf(filename_interval_string, "%s", filename_interval);
		}

		fprintf(stderr, "\n");
		fprintf(stderr, " -----  found immutable stream \"%s\" in %s  -----\n", streamID, fname);
		fprintf(stderr, "        %-20s%s\n", "filename template:", filename_template);
		fprintf(stderr, "        %-20s%s\n", "filename interval:", filename_interval_string);

		/* NB: These clobber constants must match those in the mpas_stream_manager module! */
		iclobber = 0;
		if (clobber != NULL) {
			if (strstr(clobber, "never_modify") != NULL) {
				iclobber = 0;
				fprintf(stderr, "        %-20s%s\n", "clobber mode:", "never_modify");
			}
			else if (strstr(clobber, "append") != NULL) {
				iclobber = 1;
				fprintf(stderr, "        %-20s%s\n", "clobber mode:", "append");
			}
			else if (strstr(clobber, "truncate") != NULL) {             /* Synonym for "replace_files" */
				iclobber = 2;
				fprintf(stderr, "        %-20s%s\n", "clobber mode:", "truncate");
			}
			else if (strstr(clobber, "replace_files") != NULL) {        /* Synonym for "truncate" */
				iclobber = 2;
				fprintf(stderr, "        %-20s%s\n", "clobber mode:", "replace_files");
			}
			else if (strstr(clobber, "overwrite") != NULL) {
				iclobber = 3;
				fprintf(stderr, "        %-20s%s\n", "clobber mode:", "overwrite");
			}
			else {
				iclobber = 0;
				fprintf(stderr, "        *** unrecognized clobber_mode specification; existing files will not be modified\n");
			}
		}

		/* NB: These io_type constants must match those in the mpas_stream_manager module! */
		i_iotype = 0;
		if (iotype != NULL) {
			if (strstr(iotype, "pnetcdf,cdf5") != NULL) {
				i_iotype = 1;
				fprintf(stderr, "        %-20s%s\n", "I/O type:", "Parallel-NetCDF (CDF-5, large variable support)");
			}
			else if (strstr(iotype, "pnetcdf") != NULL) {
				i_iotype = 0;
				fprintf(stderr, "        %-20s%s\n", "I/O type:", "Parallel-NetCDF");
			}
			else if (strstr(iotype, "netcdf4") != NULL) {
				i_iotype = 3;
				fprintf(stderr, "        %-20s%s\n", "I/O type:", "NetCDF-4/HDF5");
			}
			else if (strstr(iotype, "netcdf") != NULL) {
				i_iotype = 2;
				fprintf(stderr, "        %-20s%s\n", "I/O type:", "Serial NetCDF");
			}
			else {
				i_iotype = 0;
				fprintf(stderr, "        *** unrecognized io_type specification; defaulting to Parallel-NetCDF\n");
			}
		}

		/* NB: These type constants must match those in the mpas_stream_manager module! */
		if (strstr(direction, "input") != NULL && strstr(direction, "output") != NULL) {
			itype = 3;
			fprintf(stderr, "        %-20s%s\n", "direction:", "input, output");
		}
		else if (strstr(direction, "input") != NULL) {
			itype = 1;
			fprintf(stderr, "        %-20s%s\n", "direction:", "input");
		}
		else if (strstr(direction, "output") != NULL)  {
			itype = 2;
			fprintf(stderr, "        %-20s%s\n", "direction:", "output");
		}
		else  {
			itype = 4;
			fprintf(stderr, "        %-20s%s\n", "direction:", "none");
		}

		if (reference_time != NULL) {
			snprintf(ref_time_local, 256, "%s", reference_time);
		}
		else {
			snprintf(ref_time_local, 256, "initial_time");
		}
		fprintf(stderr, "        %-20s%s\n", "reference time:", ref_time_local);

		if (record_interval != NULL) {
			snprintf(rec_intv_local, 256, "%s", record_interval);
			fprintf(stderr, "        %-20s%s\n", "record interval:", rec_intv_local);
		}
		else {
			snprintf(rec_intv_local, 256, "none");
			fprintf(stderr, "        %-20s%s\n", "record interval:", "-");
		}

		if (precision != NULL && strstr(precision, "single") != NULL) {
			iprec = 4;
		}
		else if (precision != NULL && strstr(precision, "double") != NULL) {
			iprec = 8;
		}
		else {
			iprec = 0;
			if (precision != NULL)
				fprintf(stderr, "        *** unrecognized precision specification; reverting to native precision\n");
		}
		if (iprec != 0) {
			fprintf(stderr, "        %-20s%i %s\n", "real precision:", iprec, "bytes");
		}


		/* For output streams, build the directory structure where files will be written */
		if (itype == 2 || itype == 3) {
			err = build_stream_path(streamID, filename_template, mpi_comm);
			if (err != 0) {
				*status = 1;
				return;
			}
		}

		stream_mgr_create_stream_c(manager, streamID, &itype, filename_template, filename_interval_string, ref_time_local, rec_intv_local,
					&immutable, &iprec, &iclobber, &i_iotype, &err);
		if (err != 0) {
			*status = 1;
			return;
		}

		/* Possibly add an input alarm for this stream */
		if (itype == 3 || itype == 1) {
			stream_mgr_add_alarm_c(manager, streamID, "input", "start", interval_in2, &err);
			if (err != 0) {
				*status = 1;
				return;
			}
			if ( strcmp(interval_in, interval_in2) != 0 ) {
				fprintf(stderr, "        %-20s%s (%s)\n", "input alarm:", interval_in, interval_in2);
			} else {
				fprintf(stderr, "        %-20s%s\n", "input alarm:", interval_in);
			}
		}

		/* Possibly add an output alarm for this stream */
		if (itype == 3 || itype == 2) {
			stream_mgr_add_alarm_c(manager, streamID, "output", "start", interval_out2, &err);
			if (err != 0) {
				*status = 1;
				return;
			}
			if ( strcmp(interval_out, interval_out2) != 0 ) {
				fprintf(stderr, "        %-20s%s (%s)\n", "output alarm:", interval_out, interval_out2);
			} else {
				fprintf(stderr, "        %-20s%s\n", "output alarm:", interval_out);
			}
		}

		/* Possibly add packages */
		if (packagelist != NULL) {

			packages = strdup(packagelist);
			package = strsep(&packages, ";");

			stream_mgr_add_pkg_c(manager, streamID, package, &err);
			if (err != 0) {
				snprintf(msgbuf, MSGSIZE, "definition of stream \"%s\" references unrecognized package \"%s\".", streamID, package);
				fmt_warn(msgbuf);
			}
			else {
				fprintf(stderr, "        %-20s%s\n", "package:", package);
			}

			while ((package = strsep(&packages, ";")) != NULL) {
				stream_mgr_add_pkg_c(manager, streamID, package, &err);
				if (err != 0) {
					snprintf(msgbuf, MSGSIZE, "definition of stream \"%s\" references unrecognized package \"%s\".", streamID, package);
					fmt_warn(msgbuf);
				}
				else {
					fprintf(stderr, "        %-20s%s\n", "package:", package);
				}
			}

			free(packages);
		}
	}

	/* Next, handle modifications to mutable streams as well as new stream definitions */
	immutable = 0;
	for (stream_xml = ezxml_child(streams, "stream"); stream_xml; stream_xml = ezxml_next(stream_xml)) {
		streamID = ezxml_attr(stream_xml, "name");
		direction = ezxml_attr(stream_xml, "type");
		filename_template = ezxml_attr(stream_xml, "filename_template");
		filename_interval = ezxml_attr(stream_xml, "filename_interval");
		interval_in = ezxml_attr(stream_xml, "input_interval");
		interval_in2 = ezxml_attr(stream_xml, "input_interval");
		interval_out = ezxml_attr(stream_xml, "output_interval");
		interval_out2 = ezxml_attr(stream_xml, "output_interval");
		reference_time = ezxml_attr(stream_xml, "reference_time");
		record_interval = ezxml_attr(stream_xml, "record_interval");
		precision = ezxml_attr(stream_xml, "precision");
		packagelist = ezxml_attr(stream_xml, "packages");
		clobber = ezxml_attr(stream_xml, "clobber_mode");
		iotype = ezxml_attr(stream_xml, "io_type");

		/* Extract the input interval, if it refer to other streams */
		if ( interval_in ) {
			sprintf(interval_type, "input_interval");
			*status = extract_stream_interval(interval_in, interval_type, &interval_in2, streamID, streams);
			if ( *status != 0 ) {
				return;
			}
		}

		/* Extract the output interval, if it refer to other streams */
		if ( interval_out ) {
			sprintf(interval_type, "output_interval");
			*status = extract_stream_interval(interval_out, interval_type, &interval_out2, streamID, streams);
			if ( *status != 0 ) {
				return;
			}
		}

		/* Setup filename_interval correctly.
		 *
		 * If filename_interval is not explicitly set...
		 *  - Default to input interval if it is set and an interval (for input, and input;output streams)
		 *  - Default to output interval if it is set and an interval (for output streams)
		 *  - Default to none for all other cases
		 *
		 * After this check is complete, if filename_interval still has a value of NULL, it will be replaced with 'none'
		 */
		if ( filename_interval == NULL){
			/* Check for an input;output stream. Handle first as this is the most complicated case. */
			if ( strstr(direction, "input") != NULL && strstr(direction, "output") != NULL ) {

				/* If input interval is an interval (i.e. not initial_only or none) set filename_interval to the interval. */
				if ( strstr(interval_in, "initial_only") == NULL && strstr(interval_in, "none") == NULL ){
					filename_interval = interval_in2;

				/* If output interval is an interval (i.e. not initial_only or none) set filename_interval to the interval. */
				} else if ( strstr(interval_out, "initial_only") == NULL && strstr(interval_out, "none") == NULL ){
					filename_interval = interval_out2;
				}
			/* Check for an input stream. */
			} else if ( strstr(direction, "input") != NULL ) {
				if ( strstr(interval_in, "initial_only") == NULL && strstr(interval_in, "none") == NULL ){
					filename_interval = interval_in2;
				}

			/* Check for an output stream. */
			} else if ( strstr(direction, "output") != NULL ) {
				if ( strstr(interval_out, "initial_only") == NULL && strstr(interval_out, "none") == NULL ){
					filename_interval = interval_out2;
				}
			}
		} else {
			/* Handle the case where filename_interval has a value of either:
			 * output_interval -- Overwrite filename_interval with the interval provided in output_interval
			 * input_interval -- Overwrite filename_interval with the interval provided in input_interval
			 *
			 * In either of these cases, if the intervals are set to be initial_only or none, nullify filename_interval
			 * to force it's value to be none as well.
			 */
			if ( strstr(filename_interval, "input_interval") != NULL ) {
				if ( strstr(interval_in, "initial_only") == NULL && strstr(interval_in, "none") == NULL ) {
					filename_interval = interval_in2;
				} else {
					filename_interval = NULL;
				}
			} else if ( strstr(filename_interval, "output_interval") != NULL ) {
				if ( strstr(interval_out, "initial_only") == NULL && strstr(interval_out, "none") == NULL ) {
					filename_interval = interval_out2;
				} else {
					filename_interval = NULL;
				}
			}
		}

		if ( filename_interval == NULL ) {
			sprintf(filename_interval_string, "none");
		} else {
			sprintf(filename_interval_string, "%s", filename_interval);
		}

		fprintf(stderr, "\n");
		fprintf(stderr, " -----  found stream \"%s\" in %s  -----\n", streamID, fname);
		fprintf(stderr, "        %-20s%s\n", "filename template:", filename_template);
		fprintf(stderr, "        %-20s%s\n", "filename interval:", filename_interval_string);

		/* NB: These clobber constants must match those in the mpas_stream_manager module! */
		iclobber = 0;
		if (clobber != NULL) {
			if (strstr(clobber, "never_modify") != NULL) {
				iclobber = 0;
				fprintf(stderr, "        %-20s%s\n", "clobber mode:", "never_modify");
			}
			else if (strstr(clobber, "append") != NULL) {
				iclobber = 1;
				fprintf(stderr, "        %-20s%s\n", "clobber mode:", "append");
			}
			else if (strstr(clobber, "truncate") != NULL) {             /* Synonym for "replace_files" */
				iclobber = 2;
				fprintf(stderr, "        %-20s%s\n", "clobber mode:", "truncate");
			}
			else if (strstr(clobber, "replace_files") != NULL) {        /* Synonym for "truncate" */
				iclobber = 2;
				fprintf(stderr, "        %-20s%s\n", "clobber mode:", "replace_files");
			}
			else if (strstr(clobber, "overwrite") != NULL) {
				iclobber = 3;
				fprintf(stderr, "        %-20s%s\n", "clobber mode:", "overwrite");
			}
			else {
				iclobber = 0;
				fprintf(stderr, "        *** unrecognized clobber_mode specification; existing files will not be modified\n");
			}
		}

		/* NB: These io_type constants must match those in the mpas_stream_manager module! */
		i_iotype = 0;
		if (iotype != NULL) {
			if (strstr(iotype, "pnetcdf,cdf5") != NULL) {
				i_iotype = 1;
				fprintf(stderr, "        %-20s%s\n", "I/O type:", "Parallel-NetCDF (CDF-5, large variable support)");
			}
			else if (strstr(iotype, "pnetcdf") != NULL) {
				i_iotype = 0;
				fprintf(stderr, "        %-20s%s\n", "I/O type:", "Parallel-NetCDF");
			}
			else if (strstr(iotype, "netcdf4") != NULL) {
				i_iotype = 3;
				fprintf(stderr, "        %-20s%s\n", "I/O type:", "NetCDF-4/HDF5");
			}
			else if (strstr(iotype, "netcdf") != NULL) {
				i_iotype = 2;
				fprintf(stderr, "        %-20s%s\n", "I/O type:", "Serial NetCDF");
			}
			else {
				i_iotype = 0;
				fprintf(stderr, "        *** unrecognized io_type specification; defaulting to Parallel-NetCDF\n");
			}
		}

		/* NB: These type constants must match those in the mpas_stream_manager module! */
		if (strstr(direction, "input") != NULL && strstr(direction, "output") != NULL) {
			itype = 3;
			fprintf(stderr, "        %-20s%s\n", "direction:", "input, output");
		}
		else if (strstr(direction, "input") != NULL) {
			itype = 1;
			fprintf(stderr, "        %-20s%s\n", "direction:", "input");
		}
		else if (strstr(direction, "output") != NULL)  {
			itype = 2;
			fprintf(stderr, "        %-20s%s\n", "direction:", "output");
		}
		else  {
			itype = 4;
			fprintf(stderr, "        %-20s%s\n", "direction:", "none");
		}

		if (reference_time != NULL) {
			snprintf(ref_time_local, 256, "%s", reference_time);
		}
		else {
			snprintf(ref_time_local, 256, "initial_time");
		}
		fprintf(stderr, "        %-20s%s\n", "reference time:", ref_time_local);

		if (record_interval != NULL) {
			snprintf(rec_intv_local, 256, "%s", record_interval);
			fprintf(stderr, "        %-20s%s\n", "record interval:", rec_intv_local);
		}
		else {
			snprintf(rec_intv_local, 256, "none");
			fprintf(stderr, "        %-20s%s\n", "record interval:", "-");
		}

		if (precision != NULL && strstr(precision, "single") != NULL) {
			iprec = 4;
		}
		else if (precision != NULL && strstr(precision, "double") != NULL) {
			iprec = 8;
		}
		else {
			iprec = 0;
			if (precision != NULL)
				fprintf(stderr, "        *** unrecognized precision specification; reverting to native precision\n");
		}
		if (iprec != 0) {
			fprintf(stderr, "        %-20s%i %s\n", "real precision:", iprec, "bytes");
		}


		/* For output streams, build the directory structure where files will be written */
		if (itype == 2 || itype == 3) {
			err = build_stream_path(streamID, filename_template, mpi_comm);
			if (err != 0) {
				*status = 1;
				return;
			}
		}

		stream_mgr_create_stream_c(manager, streamID, &itype, filename_template, filename_interval_string, ref_time_local, rec_intv_local,
					&immutable, &iprec, &iclobber, &i_iotype, &err);
		if (err != 0) {
			*status = 1;
			return;
		}

		/* Possibly add an input alarm for this stream */
		if (itype == 3 || itype == 1) {
			stream_mgr_add_alarm_c(manager, streamID, "input", "start", interval_in2, &err);
			if (err != 0) {
				*status = 1;
				return;
			}
			if ( strcmp(interval_in, interval_in2) != 0 ) {
				fprintf(stderr, "        %-20s%s (%s)\n", "input alarm:", interval_in, interval_in2);
			} else {
				fprintf(stderr, "        %-20s%s\n", "input alarm:", interval_in);
			}
		}

		/* Possibly add an output alarm for this stream */
		if (itype == 3 || itype == 2) {
			stream_mgr_add_alarm_c(manager, streamID, "output", "start", interval_out2, &err);
			if (err != 0) {
				*status = 1;
				return;
			}
			if ( strcmp(interval_out, interval_out2) != 0 ) {
				fprintf(stderr, "        %-20s%s (%s)\n", "output alarm:", interval_out, interval_out2);
			} else {
				fprintf(stderr, "        %-20s%s\n", "output alarm:", interval_out);
			}
		}

		/* Possibly add packages */
		if (packagelist != NULL) {

			packages = strdup(packagelist);
			package = strsep(&packages, ";");

			stream_mgr_add_pkg_c(manager, streamID, package, &err);
			if (err != 0) {
				snprintf(msgbuf, MSGSIZE, "definition of stream \"%s\" references unrecognized package \"%s\".", streamID, package);
				fmt_warn(msgbuf);
			}
			else {
				fprintf(stderr, "        %-20s%s\n", "package:", package);
			}

			while ((package = strsep(&packages, ";")) != NULL) {
				stream_mgr_add_pkg_c(manager, streamID, package, &err);
				if (err != 0) {
					snprintf(msgbuf, MSGSIZE, "definition of stream \"%s\" references unrecognized package \"%s\".", streamID, package);
					fmt_warn(msgbuf);
				}
				else {
					fprintf(stderr, "        %-20s%s\n", "package:", package);
				}
			}

			free(packages);
		}

		for (varfile_xml = ezxml_child(stream_xml, "file"); varfile_xml; varfile_xml = ezxml_next(varfile_xml)) {
			varfile = ezxml_attr(varfile_xml, "name");
			packagelist = ezxml_attr(varfile_xml, "packages");

			if (packagelist != NULL)
				strncpy(packages_local, packagelist, (size_t)256);
			else
				packages_local[0] = '\0';

			/* TODO: We should probably only have one task open and read the file... */
			/* TODO: This doesn't seem like it supports var_arrays, var_structs, or streams.... */
			fd = fopen(varfile, "r");
			if (fd != NULL) {
				while (fscanf(fd, "%s", fieldname) != EOF) {
					stream_mgr_add_field_c(manager, streamID, (const char *)fieldname, packages_local, &err);
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
			packagelist = ezxml_attr(var_xml, "packages");

			if (packagelist != NULL)
				strncpy(packages_local, packagelist, (size_t)256);
			else
				packages_local[0] = '\0';

			stream_mgr_add_field_c(manager, streamID, fieldname_const, packages_local, &err);
			if (err != 0) {
				*status = 1;
				return;
			}
		}

		for (vararray_xml = ezxml_child(stream_xml, "var_array"); vararray_xml; vararray_xml = ezxml_next(vararray_xml)) {
			fieldname_const = ezxml_attr(vararray_xml, "name");
			packagelist = ezxml_attr(vararray_xml, "packages");

			if (packagelist != NULL)
				strncpy(packages_local, packagelist, (size_t)256);
			else
				packages_local[0] = '\0';

			stream_mgr_add_field_c(manager, streamID, fieldname_const, packages_local, &err);
			if (err != 0) {
				*status = 1;
				return;
			}
		}
		for (varstruct_xml = ezxml_child(stream_xml, "var_struct"); varstruct_xml; varstruct_xml = ezxml_next(varstruct_xml)) {
			structname_const = ezxml_attr(varstruct_xml, "name");
			packagelist = ezxml_attr(varstruct_xml, "packages");

			if (packagelist != NULL)
				strncpy(packages_local, packagelist, (size_t)256);
			else
				packages_local[0] = '\0';

			stream_mgr_add_pool_c(manager, streamID, structname_const, packages_local, &err);
			if (err != 0){
				*status = 1;
				return;
			}
		}

		for (substream_xml = ezxml_child(stream_xml, "stream"); substream_xml; substream_xml = ezxml_next(substream_xml)) {
			streamname_const = ezxml_attr(substream_xml, "name");
			packagelist = ezxml_attr(substream_xml, "packages");

			if (packagelist != NULL)
				strncpy(packages_local, packagelist, (size_t)256);
			else
				packages_local[0] = '\0';


			/* Immutable streams are added through the stream_mgr_add_immutable_stream_fields_c function, since
			 * they aren't defined in the XML file, and are instead defined in Registry.xml.
			 */
			stream_mgr_add_immutable_stream_fields_c(manager, streamID, streamname_const, packages_local, &err);
			if (err != 0) {
				/* If that call was successful, we DID add an immutable_stream, so do not attempt to add
				 * a mutable stream below.  Otherwise do attempt to add a mutable stream and continue.
				 */

				for(streammatch_xml = ezxml_child(streams, "stream"); streammatch_xml; streammatch_xml = ezxml_next(streammatch_xml)) {
					compstreamname_const = ezxml_attr(streammatch_xml, "name");

					if (strcmp(streamname_const, compstreamname_const) == 0) {
						for (var_xml = ezxml_child(streammatch_xml, "var"); var_xml; var_xml = ezxml_next(var_xml)) {
							fieldname_const = ezxml_attr(var_xml, "name");
							stream_mgr_add_field_c(manager, streamID, fieldname_const, packages_local, &err);
							if (err != 0) {
								*status = 1;
								return;
							}
						}


						for (vararray_xml = ezxml_child(streammatch_xml, "var_array"); vararray_xml; vararray_xml = ezxml_next(vararray_xml)) {
							fieldname_const = ezxml_attr(vararray_xml, "name");
							stream_mgr_add_field_c(manager, streamID, fieldname_const, packages_local, &err);
							if (err != 0) {
								*status = 1;
								return;
							}
						}

						for (varstruct_xml = ezxml_child(streammatch_xml, "var_struct"); varstruct_xml; varstruct_xml = ezxml_next(varstruct_xml)) {
							structname_const = ezxml_attr(varstruct_xml, "name");
							stream_mgr_add_pool_c(manager, streamID, structname_const, packages_local, &err);
							if (err != 0){
								*status = 1;
								return;
							}
						}
					}
				}
			}
		}
	}

	free(xml_buf);

	fprintf(stderr, "\n");
	fprintf(stderr, " ----- done parsing run-time I/O from %s -----\n\n", fname);
}


/*********************************************************************************
 *
 *  Function: xml_stream_get_attribute
 *
 *  Parses an XML file and searches for the stream whose name matches the 'streamname'
 *  argument; then, returns the associated attributes for that stream in
 *  the 'filename', 'ref_time', 'filename_interval', and 'io_type' arguments.
 *
 *  The fname argument provides the name of the XML file that contains the stream
 *  definitions, and mpi_comm is the Fortran MPI communicator used by MPAS.
 *
 *********************************************************************************/
void xml_stream_get_attributes(char *fname, char *streamname, int *mpi_comm, char *filename, char *ref_time, char *filename_interval, char *io_type, int *status)
{
	char *xml_buf;
	size_t bufsize;
	ezxml_t streams;
	ezxml_t stream_xml;
	const char *streamID, *filename_template, *reference_time, *c_filename_interval, *xml_iotype;
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
		reference_time = ezxml_attr(stream_xml, "reference_time");
		c_filename_interval = ezxml_attr(stream_xml, "filename_interval");
		xml_iotype = ezxml_attr(stream_xml, "io_type");

		if (strcmp(streamID, streamname) == 0) {
			found = 1;
			fprintf(stderr, "Found mesh stream with filename template %s\n", filename_template);
			sprintf(filename, "%s", filename_template);
			if ( reference_time == NULL ) {
				sprintf(ref_time, "initial_time");
			} else {
				sprintf(ref_time, "%s", reference_time);
			}

			if ( c_filename_interval == NULL ) {
				sprintf(filename_interval, "none");
			} else if ( strstr(c_filename_interval, "interval") ) {
				sprintf(filename_interval, "%s", c_filename_interval);
				c_filename_interval = ezxml_attr(stream_xml, filename_interval);

				if ( c_filename_interval == NULL ) {
					sprintf(filename_interval, "none");
				} else {
					sprintf(filename_interval, "%s", c_filename_interval);
				}
			} else {
				sprintf(filename_interval, "%s", c_filename_interval);
			}

			if ( xml_iotype == NULL ) {
				fprintf(stderr, "Using default io_type for mesh stream\n");
				sprintf(io_type, "pnetcdf");
			} else {
				if (strstr(xml_iotype, "pnetcdf,cdf5") != NULL) {
					sprintf(io_type, "%s", xml_iotype);
					fprintf(stderr, "Using io_type Parallel-NetCDF (CDF-5, large variable support) for mesh stream\n");
				}
				else if (strstr(xml_iotype, "pnetcdf") != NULL) {
					sprintf(io_type, "%s", xml_iotype);
					fprintf(stderr, "Using io_type Parallel-NetCDF for mesh stream\n");
				}
				else if (strstr(xml_iotype, "netcdf4") != NULL) {
					sprintf(io_type, "%s", xml_iotype);
					fprintf(stderr, "Using io_type NetCDF-4/HDF5 for mesh stream\n");
				}
				else if (strstr(xml_iotype, "netcdf") != NULL) {
					sprintf(io_type, "%s", xml_iotype);
					fprintf(stderr, "Using io_type Serial NetCDF for mesh stream\n");
				}
				else {
					sprintf(io_type, "pnetcdf");
					fprintf(stderr, "*** unrecognized io_type specification for mesh stream; defaulting to Parallel-NetCDF\n");
				}
			}
			break;
		}
	}
	if (found == 0) {
		*status = 1;
		return;
	}
}
