// Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
// and the University Corporation for Atmospheric Research (UCAR).
//
// Unless noted otherwise source code is licensed under the BSD license.
// Additional copyright and license information can be found in the LICENSE file
// distributed with this code, or at http://mpas-dev.github.com/license.html
//
#define INTEGER 0
#define REAL 1
#define LOGICAL 2
#define CHARACTER 3

#define PERSISTENT 0
#define SCRATCH    1

#define INPUT0    0x00000001
#define RESTART0  0x00000002
#define OUTPUT0   0x00000004
#define SFC0      0x00000008

#define NEW_NAMELIST(X) X = (struct namelist *)malloc(sizeof(struct namelist)); X->next = NULL;
#define NEW_DIMENSION(X) X = (struct dimension *)malloc(sizeof(struct dimension)); X->next = NULL;
#define NEW_DIMENSION_LIST(X) X = (struct dimension_list *)malloc(sizeof(struct dimension_list)); X->dim = NULL; X->prev = NULL; X->next = NULL;
#define NEW_VARIABLE(X) X = (struct variable *)malloc(sizeof(struct variable)); X->dimlist = NULL; X->next = NULL;
#define NEW_VARIABLE_LIST(X) X = (struct variable_list *)malloc(sizeof(struct variable_list)); X->var = NULL; X->prev = NULL; X->next = NULL;
#define NEW_GROUP_LIST(X) X = (struct group_list *)malloc(sizeof(struct group_list)); X->vlist = NULL; X->next = NULL;

union default_val {
   int ival;
   float rval;
   int lval;
   char cval[32];
};

struct namelist {
   char name[1024];
   char record[1024];
   int vtype;
   union default_val defval;
   struct namelist * next;
};

struct dimension {
   char name_in_file[1024];
   char name_in_code[1024];
   int constant_value;
   int namelist_defined;
   struct dimension * next;
};

struct dimension_list {
   struct dimension * dim;
   struct dimension_list * prev;
   struct dimension_list * next;
};

struct variable_list {
   struct variable * var;
   struct variable_list * prev;
   struct variable_list * next;
};

struct group_list {
   char name[1024];
   struct variable_list * vlist; 
   struct group_list * next; 
   int ntime_levs;
};

struct variable {
   char name_in_file[1024];
   char name_in_code[1024];
   char struct_group[1024];
   char super_array[1024];
   char array_class[1024];
   char default_value[1024];
   int persistence;
   int vtype;
   int ndims;
   int timedim;
   int iostreams;
   int decomposed;
   struct dimension_list * dimlist;
   struct variable * next;
};
