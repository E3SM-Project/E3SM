// Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
// and the University Corporation for Atmospheric Research (UCAR).
//
// Unless noted otherwise source code is licensed under the BSD license.
// Additional copyright and license information can be found in the LICENSE file
// distributed with this code, or at http://mpas-dev.github.com/license.html
//

#include <string.h>
#include "ezxml.h"

int is_structure_writable(ezxml_t structure, int pairs, char **keys, char **values){/*{{{*/
	const char *control_value;

	int i, value_found;

	char *string;
	char *tofree;
	char *token;

	value_found = -1;

	for (i = 0; i < pairs; i++){
		control_value = ezxml_attr(structure, keys[i]);

		if ( control_value != NULL ){
			/* Set found to false, since the key was found. */
			value_found = 0;
			string = strdup(control_value);
			tofree = string;

			while( ( token = strsep(&string, ";") ) != NULL ){
				if ( strcmp(token, values[i]) == 0 ){
					value_found = 1;
				}
			}

		}
	}

	return value_found;
}/*}}}*/
