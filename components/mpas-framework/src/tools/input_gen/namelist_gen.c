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

int generate_namelist(ezxml_t registry, FILE* fd, int pairs, char **keys, char **values);
char* get_option_value(ezxml_t nml_option, int pairs, char **values);

int main(int argc, char ** argv)/*{{{*/
{
	FILE *regfile, *namelistfile;
	int i, err, pairs, pairs_allc;
	char **keys, **values;
	char *string, *tofree, *token;

	if (! (argc >= 3) ) {
		fprintf(stderr, "\nError: Not enough input arguments...\n\n");
		fprintf(stderr, "\nUSAGE: ./namelist_gen registry_file namelist_file_output [value_attribute] [key1=value1] [key2=value2]...\n\n");
		return 1;
	}

	if (!(regfile = fopen(argv[1], "r"))) {
		fprintf(stderr,"\nError: Could not open file %s for reading.\n\n", argv[1]);
		return 1;
	}
	if (!(namelistfile = fopen(argv[2], "w+"))) {
		fprintf(stderr, "\nError: Could not create file %s for writing.\n\n", argv[2]);
		return 1;
	}

	pairs = argc-3;
	if ( pairs > 0 ) {
		keys = malloc(sizeof(char*) * pairs);
		values = malloc(sizeof(char*) * pairs);
	} else {
		pairs = 0;
	}

	for(i = 0; i < pairs; i++){
		keys[i] = malloc(sizeof(char)*1024);
		values[i] = malloc(sizeof(char)*1024);

		string = strdup(argv[3+i]);
		tofree = string;

		token = strsep(&string, "=");
		sprintf(keys[i], "%s", token);

		token = strsep(&string, "=");
		sprintf(values[i], "%s", token);

		free(tofree);
	}

	ezxml_t registry = ezxml_parse_fp(regfile);

	err = generate_namelist(registry, namelistfile, pairs, keys, values);

	if ( pairs > 0 ) {
		free(keys);
		free(values);
	}

	fclose(regfile);
	fclose(namelistfile);

	return 0;
}/*}}}*/

int generate_namelist(ezxml_t registry, FILE* fd, int pairs, char **keys, char **values){/*{{{*/
	ezxml_t nmlrecs_xml, nmlopt_xml;

	const char *recname;
	const char *optname, *opttype;
	char *optval;

	int write_opt, write_record, record_written;
	int i;

	for (nmlrecs_xml = ezxml_child(registry, "nml_record"); nmlrecs_xml; nmlrecs_xml = nmlrecs_xml->next){
		recname = ezxml_attr(nmlrecs_xml, "name");

		record_written = 0;
		write_record = is_structure_writable(nmlrecs_xml, pairs, keys, values);

		/* If record doesn't have a key, record is writeable */
		if ( write_record == -1 ) write_record = 1;

		if ( write_record == 1 ) {
			fprintf(fd, "&%s\n", recname);
			record_written = 1;
		}

		for (nmlopt_xml = ezxml_child(nmlrecs_xml, "nml_option"); nmlopt_xml; nmlopt_xml = nmlopt_xml->next){
			optname = ezxml_attr(nmlopt_xml, "name");
			opttype = ezxml_attr(nmlopt_xml, "type");

			write_opt = is_structure_writable(nmlopt_xml, pairs, keys, values);

			if ( write_opt == 1 || ( write_opt == -1 && write_record == 1) ){
				if ( !record_written ) {
					fprintf(fd, "&%s\n", recname);
					record_written = 1;
				}

				optval = get_option_value(nmlopt_xml, pairs, values);

				if ( strcmp(opttype, "character") == 0){
					fprintf(fd, "    %s = '%s'\n", optname, optval);
				} else {
					fprintf(fd, "    %s = %s\n", optname, optval);
				}

				free(optval);
			}
		}

		if ( record_written ) {
			fprintf(fd, "/\n");
		}
	}

	return 0;
}/*}}}*/

char* get_option_value(ezxml_t nml_option, int pairs, char **values){/*{{{*/
	int i;
	char value_name[1024];
	const char *test_val;
	char *optvalue;

	for (i = 0; i < pairs; i++){
		sprintf(value_name, "%s_value", values[i]);

		test_val = ezxml_attr(nml_option, value_name);

		if ( test_val != NULL ) {
			optvalue = strdup(test_val);
			return optvalue;
		}
	}

	test_val = ezxml_attr(nml_option, "default_value");
	optvalue = strdup(test_val);
	return optvalue;
}/*}}}*/
