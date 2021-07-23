// Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
// and the University Corporation for Atmospheric Research (UCAR).
//
// Unless noted otherwise source code is licensed under the BSD license.
// Additional copyright and license information can be found in the LICENSE file
// distributed with this code, or at http://mpas-dev.github.com/license.html
//

#include <stdlib.h>
#include <time.h>

/* Use the following interface in Fortran for seed_random()

 interface
   subroutine seed_random() bind(c)
   end subroutine seed_random
 end interface

*/
void seed_random() {
	srand(time(NULL));
}

/* Use the following interface in Fortran for gen_random()

  interface
    subroutine gen_random(len, id) bind(c)
       use iso_c_binding, only : c_int, c_char
       integer (c_int), intent(in), value :: len
       character (c_char), dimension(*), intent(inout) :: id
    end subroutine gen_random
  end interface

*/
void gen_random(int len, char * id) {/*{{{*/
	int i;
	int r;
	static const char alphanum[] =
		"0123456789"
		"abcdefghijklmnopqrstuvwxyz";

	for (i = 0; i < len; ++i) {
		r = rand();
		id[i] = alphanum[r % (sizeof(alphanum) - 1)]; 
	}    

}/*}}}*/
