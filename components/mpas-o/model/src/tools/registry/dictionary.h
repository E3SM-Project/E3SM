// Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
// and the University Corporation for Atmospheric Research (UCAR).
//
// Unless noted otherwise source code is licensed under the BSD license.
// Additional copyright and license information can be found in the LICENSE file
// distributed with this code, or at http://mpas-dev.github.com/license.html
//

#define TABLESIZE 271

struct dnode {
	char key[1024];
	struct dnode * next;
};

struct dtable {
	int size;
	struct dnode * table[TABLESIZE];
};

void dict_alloc(struct dtable **);
void dict_insert(struct dtable *, char *);
void dict_remove(struct dtable *, char *);
int dict_search(struct dtable *, char *);
int dict_size(struct dtable *);
void dict_free(struct dtable **);
