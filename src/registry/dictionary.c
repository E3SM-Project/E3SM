// Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
// and the University Corporation for Atmospheric Research (UCAR).
//
// Unless noted otherwise source code is licensed under the BSD license.
// Additional copyright and license information can be found in the LICENSE file
// distributed with this code, or at http://mpas-dev.github.com/license.html
//
#include <stdlib.h>
#include <string.h>
#include "dictionary.h"

int hashstring(char *);

void dict_alloc(struct dtable ** dict)
{
	int i;

	*dict = (struct dtable *)malloc(sizeof(struct dtable));

	for(i=0; i<TABLESIZE; i++)
		(*dict)->table[i] = NULL;

	(*dict)->size = 0;
}


void dict_insert(struct dtable * dict, char * word)
{
	int hval;
	struct dnode * dptr;

	hval = hashstring(word) % TABLESIZE;

	dptr = (struct dnode *)malloc(sizeof(struct dnode));
	strncpy(dptr->key, word, 1024);
	dptr->next = dict->table[hval];
	dict->table[hval] = dptr;

	dict->size++;
}


void dict_remove(struct dtable * dict, char * word)
{
	int hval;
	struct dnode * dptr_prev;
	struct dnode * dptr;

	hval = hashstring(word) % TABLESIZE;

	dptr_prev = 0;
	dptr = dict->table[hval];

	while (dptr && strncmp(dptr->key, word, 1024) != 0) {
		dptr_prev = dptr;
		dptr = dptr->next;
	}

	if (dptr) {
		if (dptr_prev)
			dptr_prev->next = dptr->next;
		else
			dict->table[hval] = dict->table[hval]->next;
		free(dptr);
		dict->size--;
	}
}


int dict_search(struct dtable * dict, char * word)
{
	int hval;
	struct dnode * dptr;

	hval = hashstring(word) % TABLESIZE;

	dptr = dict->table[hval];
	while (dptr && strncmp(dptr->key, word, 1024) != 0)
		dptr = dptr->next;

	if (!dptr) return 0;

	return 1;
}


int dict_size(struct dtable * dict)
{
	return dict->size;
}


void dict_free(struct dtable ** dict)
{
	int i;
	struct dnode * dptr;

	for(i=0; i<TABLESIZE; i++) {
		while ((*dict)->table[i]) {
			dptr = (*dict)->table[i];
			(*dict)->table[i] = (*dict)->table[i]->next;
			free(dptr);
		}
	}

	free(*dict);
}


int hashstring(char * word)
{
	int i;
	int hval;

	hval = 0;

	for(i=0; i<1024 && word[i] != '\0'; i++) {
		hval = hval + (int)word[i];
	}

	return hval;
}
