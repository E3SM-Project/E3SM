// Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
// and the University Corporation for Atmospheric Research (UCAR).
//
// Unless noted otherwise source code is licensed under the BSD license.
// Additional copyright and license information can be found in the LICENSE file
// distributed with this code, or at http://mpas-dev.github.com/license.html
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ezxml.h"
#include "test_functions.h"

#define SINGLE 0
#define SEPARATE 1

#define ORDER_LISTED 0
#define ORDER_MUTABLE 1

int generate_streams(ezxml_t registry, FILE* fd, char *stream_file_prefix, int order, int pairs, char **keys, char **values);
void write_stream_header(ezxml_t stream_xml, FILE *fd);

int main(int argc, char ** argv)/*{{{*/
{
	FILE * regfile, * streamsfile;
	int i, err, pairs;
	char *stream_file_prefix;
	char **keys, **values;
	char *string, *tofree, *token;
	int order;

	if (! (argc >= 5) ) {
		fprintf(stderr, "\nError: Not enough input arguments...\n\n");
		fprintf(stderr, "\nUSAGE: ./streams_gen registry_file streams_file_output stream_file_prefix stream_order [key1=value1] [key2=value2]...\n\n");
		return 1;
	}

	if (!(regfile = fopen(argv[1], "r"))) {
		fprintf(stderr,"\nError: Could not open file %s for reading.\n\n", argv[1]);
		return 1;
	}
	if (!(streamsfile = fopen(argv[2], "w+"))) {
		fprintf(stderr, "\nError: Could not create file %s for writing.\n\n", argv[2]);
		return 1;
	}

	stream_file_prefix=strdup(argv[3]);

	if (strcmp(argv[4], "listed") == 0 ){
		order = ORDER_LISTED;
	} else if (strcmp(argv[4], "mutable") == 0){
		order = ORDER_MUTABLE;
	} else {
		fprintf(stderr, "\nError: Option %s is not valid for stream order.\n", argv[4]);
		fprintf(stderr, "\tOrder values are:\n");
		fprintf(stderr, "\t\tlisted -- Defines streams in the listed order in registry.\n");
		fprintf(stderr, "\t\tmutable -- Defines immutable streams first, followed by mutable streams.\n");
		return 1;
	}


	pairs = argc-5;
	if ( pairs > 0 ) {
		keys = malloc(sizeof(char*) * pairs);
		values = malloc(sizeof(char*) * pairs);
	} else {
		pairs = 0;
	}

	for(i = 0; i < pairs; i++){
		keys[i] = malloc(sizeof(char)*1024);
		values[i] = malloc(sizeof(char)*1024);

		string = strdup(argv[5+i]);
		tofree = string;

		token = strsep(&string, "=");
		sprintf(keys[i], "%s", token);

		token = strsep(&string, "=");
		sprintf(values[i], "%s", token);

		free(tofree);
	}

	ezxml_t registry = ezxml_parse_fp(regfile);

	err = generate_streams(registry, streamsfile, stream_file_prefix, order, pairs, keys, values);

	if ( pairs > 0 ) {
		free(keys);
		free(values);
	}
	free(stream_file_prefix);

	fclose(regfile);
	fclose(streamsfile);

	return 0;
}/*}}}*/

int generate_streams(ezxml_t registry, FILE* fd, char *stream_file_prefix, int order, int pairs, char **keys, char **values){/*{{{*/
	ezxml_t streams_xml, stream_xml;
	ezxml_t member_xml;

	const char *name, *type, *immutable, *filename_template, *filename_interval, *packages, *record_interval;
	const char *varpackages;
	const char *reference_time, *clobber_mode, *useMissingValMask, *precision, *input_interval, *output_interval;

	const char *runtime, *subname;

	char filename[256];

	int write_stream, write_member, stream_written, filetype;

	FILE *fd2;

	fprintf(fd, "<streams>\n");

	if ( order == ORDER_MUTABLE ) {
		for (streams_xml = ezxml_child(registry, "streams"); streams_xml; streams_xml = streams_xml->next){

			// First pass, only pull out immutable streams. Put these at the top of the file. Since they *have* to be defined.
			for (stream_xml = ezxml_child(streams_xml, "stream"); stream_xml; stream_xml = stream_xml->next){
				immutable = ezxml_attr(stream_xml, "immutable");

				write_stream = is_structure_writable(stream_xml, pairs, keys, values);

				if ( write_stream == -1 ) write_stream = 1; // If key is missing, make stream writable

				if ( immutable != NULL && strcmp(immutable, "true") == 0) {
					if ( write_stream == 1 ) {
						write_stream_header(stream_xml, fd);
					}
				}
			}
		}
	}

	for (streams_xml = ezxml_child(registry, "streams"); streams_xml; streams_xml = streams_xml->next){
		/* Second pass is intended to pull out mutable streams when order is mutable, but pulls out all streams when order is listed. */
		for (stream_xml = ezxml_child(streams_xml, "stream"); stream_xml; stream_xml = stream_xml->next){
			stream_written = 0;

			name = ezxml_attr(stream_xml, "name");
			runtime = ezxml_attr(stream_xml, "runtime_format");

			immutable = ezxml_attr(stream_xml, "immutable");

			write_stream = is_structure_writable(stream_xml, pairs, keys, values);

			if ( write_stream == -1 ) write_stream = 1; // If key is missing, make stream writable

			/* If order was listed, then write immutable streams where they are defined. */
			if ( order == ORDER_LISTED ) {
				if ( immutable != NULL && strcmp(immutable, "true") == 0){
					if ( write_stream == 1) {
						write_stream_header(stream_xml, fd);
					}
				}
			}

			if ( immutable == NULL || (immutable != NULL && strcmp(immutable, "false") == 0) ) {

				if ( runtime == NULL ) {
					fprintf(stderr, " Stream %s requires the runtime_format attribute.\n", name);
				}

				if ( strcmp(runtime, "single_file") == 0 ){
					filetype = SINGLE;
				} else if ( strcmp(runtime, "separate_file") == 0){
					filetype = SEPARATE;
					sprintf(filename, "%s%s", stream_file_prefix, name);
				} else {
					fprintf(stderr, " Stream %s has attribute runtime_format set to an invalid option: %s\n", name, runtime);
					return 1;
				}

				write_stream = is_structure_writable(stream_xml, pairs, keys, values);

				// If stream doesn't contain the key, make it writable anyway.
				if ( write_stream == -1 ) write_stream = 1;

				if ( write_stream ){
					write_stream_header(stream_xml, fd);
					stream_written = 1;
					if ( filetype == SEPARATE ) {
						fprintf(fd, "\t<file name=\"%s\"/>\n", filename);
						fd2 = fopen(filename, "w+");
					}
				}

				// Extract all streams from within the stream.
				for (member_xml = ezxml_child(stream_xml, "stream"); member_xml; member_xml = member_xml->next){
					write_member = is_structure_writable(member_xml, pairs, keys, values);

					if ( write_member == 1 || (write_member == -1 && write_stream) ) {
						if ( !stream_written ) {
							write_stream_header(stream_xml, fd);
							stream_written = 1;
							if ( filetype == SEPARATE ) {
								fprintf(fd, "\t<file name=\"%s\"/>\n", filename);
								fd2 = fopen(filename, "w+");
							}
						}

						subname = ezxml_attr(member_xml, "name");
						varpackages = ezxml_attr(member_xml, "packages");

						if (varpackages == NULL) {
							fprintf(fd, "\t<stream name=\"%s\"/>\n", subname);
						} else {
							fprintf(fd, "\t<stream name=\"%s\" packages=\"%s\"/>\n", subname, varpackages);
						}
					}
				}

				// Extract all var_structs from within the stream.
				for (member_xml = ezxml_child(stream_xml, "var_struct"); member_xml; member_xml = member_xml->next){
					write_member = is_structure_writable(member_xml, pairs, keys, values);

					if ( write_member == 1 || (write_member == -1 && write_stream) ) {
						if ( !stream_written ) {
							write_stream_header(stream_xml, fd);
							stream_written = 1;
							if ( filetype == SEPARATE ) {
								fprintf(fd, "\t<file name=\"%s\"/>\n", filename);
								fd2 = fopen(filename, "w+");
							}
						}

						subname = ezxml_attr(member_xml, "name");
						varpackages = ezxml_attr(member_xml, "packages");

						if (varpackages == NULL) {
							fprintf(fd, "\t<var_struct name=\"%s\"/>\n", subname);
						} else {
							fprintf(fd, "\t<var_struct name=\"%s\" packages=\"%s\"/>\n", subname, varpackages);
						}
					}
				}

				// Extract all var_arrays from within the stream.
				for (member_xml = ezxml_child(stream_xml, "var_array"); member_xml; member_xml = member_xml->next){
					write_member = is_structure_writable(member_xml, pairs, keys, values);

					if ( write_member == 1 || (write_member == -1 && write_stream) ) {
						if ( !stream_written ) {
							write_stream_header(stream_xml, fd);
							stream_written = 1;
							if ( filetype == SEPARATE ) {
								fprintf(fd, "\t<file name=\"%s\"/>\n", filename);
								fd2 = fopen(filename, "w+");
							}
						}

						subname = ezxml_attr(member_xml, "name");
						varpackages = ezxml_attr(member_xml, "packages");

						if ( filetype == SINGLE ) {
							if (varpackages == NULL) {
								fprintf(fd, "\t<var_array name=\"%s\"/>\n", subname);
							} else {
								fprintf(fd, "\t<var_array name=\"%s\" packages=\"%s\"/>\n", subname, varpackages);
							}
						} else if ( filetype == SEPARATE ) {
							fprintf(fd2, "%s\n", subname);
							if (varpackages != NULL) {
								fprintf(stderr, "Warning: Unable to add packages \"%s\" to var_array \"%s\" in stream \"%s\" because \"runtime_format=separate_file\".\n", varpackages, subname, name);
							}
						}
					}
				}

				// Extract all vars from within the stream.
				for (member_xml = ezxml_child(stream_xml, "var"); member_xml; member_xml = member_xml->next){
					write_member = is_structure_writable(member_xml, pairs, keys, values);

					if ( write_member == 1 || (write_member == -1 && write_stream) ) {
						if ( !stream_written ) {
							write_stream_header(stream_xml, fd);
							stream_written = 1;
							if ( filetype == SEPARATE ) {
								fprintf(fd, "\t<file name=\"%s\"/>\n", filename);
								fd2 = fopen(filename, "w+");
							}
						}
						subname = ezxml_attr(member_xml, "name");
						varpackages = ezxml_attr(member_xml, "packages");

						if ( filetype == SINGLE ) {
							if (varpackages == NULL) {
								fprintf(fd, "\t<var name=\"%s\"/>\n", subname);
							} else {
								fprintf(fd, "\t<var name=\"%s\" packages=\"%s\"/>\n", subname, varpackages);
							}
						} else if ( filetype == SEPARATE ) {
							fprintf(fd2, "%s\n", subname);
							if (varpackages != NULL) {
								fprintf(stderr, "Warning: Unable to add packages \"%s\" to var \"%s\" in stream \"%s\" because \"runtime_format=separate_file\".\n", varpackages, subname, name);
							}
						}
					}
				}

				if ( stream_written == 1 ) {
					fprintf(fd, "</stream>\n\n");
					if ( filetype == SEPARATE ) {
						fclose(fd2);
					}
				}
			}
		}
	}
	fprintf(fd, "</streams>\n");

	return 0;
}/*}}}*/

void write_stream_header(ezxml_t stream_xml, FILE *fd){/*{{{*/
	const char *name, *type, *immutable, *filename_template, *filename_interval, *packages, *record_interval;
	const char *reference_time, *clobber_mode, *useMissingValMask, *precision, *input_interval, *output_interval, *io_type;

	char spacing[1024];

	name = ezxml_attr(stream_xml, "name");
	type = ezxml_attr(stream_xml, "type");
	filename_template = ezxml_attr(stream_xml, "filename_template");
	filename_interval = ezxml_attr(stream_xml, "filename_interval");
	reference_time = ezxml_attr(stream_xml, "reference_time");
	clobber_mode = ezxml_attr(stream_xml, "clobber_mode");
	useMissingValMask = ezxml_attr(stream_xml, "useMissingValMask");
	input_interval = ezxml_attr(stream_xml, "input_interval");
	output_interval = ezxml_attr(stream_xml, "output_interval");
	record_interval = ezxml_attr(stream_xml, "record_interval");
	precision = ezxml_attr(stream_xml, "precision");
	io_type = ezxml_attr(stream_xml, "io_type");
	immutable = ezxml_attr(stream_xml, "immutable");
	packages = ezxml_attr(stream_xml, "packages");

	if (immutable != NULL && strcmp(immutable, "true") == 0){
		sprintf(spacing, "                  ");
		fprintf(fd, "<immutable_stream name=\"%s\"\n", name);
	} else {
		sprintf(spacing, "        ");
		fprintf(fd, "<stream name=\"%s\"\n", name);
	}

	fprintf(fd, "%stype=\"%s\"", spacing, type);

	fprintf(fd, "\n%sfilename_template=\"%s\"", spacing, filename_template);
	if ( filename_interval != NULL ){
		fprintf(fd, "\n%sfilename_interval=\"%s\"", spacing, filename_interval);
	}
	if ( reference_time != NULL ){
		fprintf(fd, "\n%sreference_time=\"%s\"", spacing, reference_time);
	}
	if ( record_interval != NULL ){
		fprintf(fd, "\n%srecord_interval=\"%s\"", spacing, record_interval);
	}
	if ( clobber_mode != NULL ){
		fprintf(fd, "\n%sclobber_mode=\"%s\"", spacing, clobber_mode);
	}
	if ( useMissingValMask != NULL ){
		fprintf(fd, "\n%suseMissingValMask=\"%s\"", spacing, useMissingValMask);
	}
	if ( precision != NULL ){
		fprintf(fd, "\n%sprecision=\"%s\"", spacing, precision);
	}
	if ( io_type != NULL ){
		fprintf(fd, "\n%sio_type=\"%s\"", spacing, io_type);
	}
	if ( packages != NULL ){
		fprintf(fd, "\n%spackages=\"%s\"", spacing, packages);
	}
	if ( input_interval != NULL ) {
		fprintf(fd, "\n%sinput_interval=\"%s\"", spacing, input_interval);
	}
	if ( output_interval != NULL ) {
		fprintf(fd, "\n%soutput_interval=\"%s\"", spacing, output_interval);
	}

	if (immutable != NULL && strcmp(immutable, "true") == 0){
		fprintf(fd, " />\n\n");
	} else {
		fprintf(fd, " >\n\n");
	}




}/*}}}*/

