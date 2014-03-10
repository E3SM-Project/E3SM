// Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
// and the University Corporation for Atmospheric Research (UCAR).
//
// Unless noted otherwise source code is licensed under the BSD license.
// Additional copyright and license information can be found in the LICENSE file
// distributed with this code, or at http://mpas-dev.github.com/license.html
//

#include <stdlib.h>
#include <time.h>

#ifdef UNDERSCORE
#define gen_random gen_random_
#else
#ifdef DOUBLEUNDERSCORE
#define gen_random gen_random__
#endif
#endif

void gen_random(int * len, char * id) {/*{{{*/
	int i;
	int r;
	static const char alphanum[] =
		"0123456789"
		"abcdefghijklmnopqrstuvwxyz";

	srand(time(NULL));

	for (i = 0; i < *len; ++i) {
		r = rand();
		id[i] = alphanum[r % (sizeof(alphanum) - 1)]; 
	}    

}/*}}}*/
