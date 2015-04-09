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
#include "ezxml/ezxml.h"
#include "registry_types.h"

int is_derived_dim(char * d)/*{{{*/
{
	if (strchr(d, (int)'+')) return 1;
	if (strchr(d, (int)'-')) return 1;

	return 0;
}/*}}}*/


char * new_dimension_name(char * old_name){/*{{{*/
	int i, j;
	int len, new_len;
	int new_string;
	int added_sections;
	int symbols;
	char * new_name;

	i = 0;
	new_string = 0;
	added_sections = 0;
	len = 0;
	symbols = 0;
	while (old_name[i] != '\0') {
		if (new_string == 0 && ((old_name[i] >= 'a' && old_name[i] <= 'z') || (old_name[i] >= 'A' && old_name[i] <= 'Z'))){
			new_string = 1;
			added_sections++;
		} else if (old_name[i] == '+' || old_name[i] == '-' || old_name[i] == '*' || old_name[i] == '/'){
			new_string = 0;
			symbols++;
		}
		len++;
		i++;
	}

	new_len = len + 7 + (added_sections-1) * 8 + symbols * 3 + 1;
	new_name = malloc(sizeof(char)*new_len);
	i = 0;
	j = 0;
	added_sections = 0;
	while (old_name[i] != '\0') {
		if (new_string == 0 && ((old_name[i] >= 'a' && old_name[i] <= 'z') || (old_name[i] >= 'A' && old_name[i] <= 'Z'))){
			new_string = 1;
			if (added_sections == 0){
				new_name[j] = 'm';
				new_name[j+1] = 'e';
				new_name[j+2] = 's';
				new_name[j+3] = 'h';
				new_name[j+4] = ' ';
				new_name[j+5] = '%';
				new_name[j+6] = ' ';

				j += 7;
			} else {
				new_name[j] = ' ';
				new_name[j+1] = 'm';
				new_name[j+2] = 'e';
				new_name[j+3] = 's';
				new_name[j+4] = 'h';
				new_name[j+5] = ' ';
				new_name[j+6] = '%';
				new_name[j+7] = ' ';

				j += 8;
			}

			added_sections++;
		} 
		if (old_name[i] == '+' || old_name[i] == '-' || old_name[i] == '*' || old_name[i] == '/'){
			new_string = 0;
			new_name[j] = ' ';

			j++;
		}

		new_name[j] = old_name[i];
		j++;

		if (old_name[i] == '+' || old_name[i] == '-' || old_name[i] == '*' || old_name[i] == '/'){
			new_string = 0;
			new_name[j] = ' ';

			j++;
		}

		i++;
	}

	new_name[j] = '\0';

	return new_name;
}/*}}}*/


void split_derived_dim_string(char * dim, char ** p1, char ** p2)/*{{{*/
{
	char * cp, * cm, * c;
	int n;

	cp = strchr(dim, (int)'+');
	cm = strchr(dim, (int)'-');
	if (!cp) 
		c = cm;
	else if (!cm) 
		c = cp;
	else if (cm < cp) 
		c = cm;
	else 
		c = cp;

	n = c - dim;
	*p1 = (char *)malloc(n*sizeof(char));
	snprintf(*p1, n, "%s", dim+1);

	*p2 = (char *)malloc((strlen(dim)-n+1)*sizeof(char));
	sprintf(*p2, "%s", dim+n);
}/*}}}*/


int is_integer_constant(char * c) {/*{{{*/
	int i;

	i = 0;
	while (c[i] != '\0') {
		if (c[i] < '0' || c[i] > '9') return -1;
		i++;
	}

	return atoi(c);
}/*}}}*/

char * check_packages(ezxml_t registry, char * packages){/*{{{*/
	ezxml_t packages_xml, package_xml;

	const char *packagename;

	char *string, *tofree, *token;
	char *failed;
	int missing_package;

	string = strdup(packages);
	tofree = string;
	failed = NULL;

	while( (token = strsep(&string, ";")) != NULL) {
		missing_package = 1;
		for (packages_xml = ezxml_child(registry, "packages"); packages_xml; packages_xml = packages_xml->next){
			for (package_xml = ezxml_child(packages_xml, "package"); package_xml; package_xml = package_xml->next){
				packagename = ezxml_attr(package_xml, "name");

				if(strcasecmp(packagename, token) == 0){
					missing_package = 0;
				}
			}
		}

		if (missing_package) {
			failed = strdup(token);
			free(tofree);
			return failed;
		}
	}
	free(tofree);
	return failed;
}/*}}}*/

char * check_dimensions(ezxml_t registry, char * dims){/*{{{*/
	ezxml_t dims_xml, dim_xml;

	const char *dimname;

	char *string, *tofree, *token;
	int missing_dim;

	string = strdup(dims);
	tofree = string;

	while( (token = strsep(&string, " ")) != NULL) {
		if (strcasecmp(token, "Time") != 0){
			missing_dim = 1;
			for (dims_xml = ezxml_child(registry, "dims"); dims_xml; dims_xml = dims_xml->next){
				for (dim_xml = ezxml_child(dims_xml, "dim"); dim_xml; dim_xml = dim_xml->next){
					dimname = ezxml_attr(dim_xml, "name");

					if(strcasecmp(dimname, token) == 0){
						missing_dim = 0;
					}
				}
			}

			if (missing_dim) {
				free(tofree);
				return token;
			}
		}
	}
	free(tofree);
	return NULL;
}/*}}}*/


char * check_streams(ezxml_t registry, char * streams)
{
	ezxml_t streams_xml, stream_xml;

	const char *streamname;

	char *string, *tofree, *token;
	char *failed;
	int missing_stream;

	string = strdup(streams);
	tofree = string;
	failed = NULL;

	while( (token = strsep(&string, ";")) != NULL) {
		missing_stream = 1;
		for (streams_xml = ezxml_child(registry, "streams"); streams_xml; streams_xml = streams_xml->next) {
			for (stream_xml = ezxml_child(streams_xml, "stream"); stream_xml; stream_xml = stream_xml->next) {
				streamname = ezxml_attr(stream_xml, "name");

				if(strcasecmp(streamname, token) == 0) {    /* TODO: Not portable? */
					missing_stream = 0;
				}
			}
		}

		if (missing_stream) {
			failed = strdup(token);
			free(tofree);
			return failed;
		}
	}
	free(tofree);
	return failed;
}


int check_persistence(const char * persistence){/*{{{*/
	if(persistence){
		if(strncmp(persistence, "persistent", 1024) == 0){
			return PERSISTENT;
		} else if(strncmp(persistence, "scratch", 1024) == 0){
			return SCRATCH;
		} else {
			fprintf(stderr, "ERROR: In check_persistence. Persistence not equal to persistent or scratch.\n");
			return -1;
		}
	} else {
		return PERSISTENT;
	}
}/*}}}*/
