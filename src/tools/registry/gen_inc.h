// Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
// and the University Corporation for Atmospheric Research (UCAR).
//
// Unless noted otherwise source code is licensed under the BSD license.
// Additional copyright and license information can be found in the LICENSE file
// distributed with this code, or at http://mpas-dev.github.com/license.html
//
//

#include "ezxml.h"

void write_model_variables(ezxml_t registry);
int write_field_pointers(FILE* fd);
int write_field_pointer_arrays(FILE* fd);
int set_pointer_name(int type, int ndims, char *pointer_name);
int add_package_to_list(const char * package, const char * package_list);
int build_struct_package_lists(ezxml_t currentPosition, char * out_packages);
int get_dimension_information(ezxml_t registry, const char *test_dimname, int *has_time, int *decomp);
int build_dimension_information(ezxml_t registry, ezxml_t var, int *ndims, int *has_time, int *decomp);
int get_field_information(const char *vartype, const char *varval, char *default_value, const char *varmissval, char *missing_value, int *type);
int write_set_field_pointer(FILE *fd, const char *spacing, const char *iterator_name, const char *pool_name);
void write_default_namelist(ezxml_t registry);
int parse_packages_from_registry(ezxml_t registry);
int parse_namelist_records_from_registry(ezxml_t registry);
int parse_dimensions_from_registry(ezxml_t registry);
int parse_var_array(FILE *fd, ezxml_t registry, ezxml_t superStruct, ezxml_t varArray, const char * corename);
int parse_var(FILE *fd, ezxml_t registry, ezxml_t superStruct, ezxml_t currentVar, const char * corename);
int parse_struct(FILE *fd, ezxml_t registry, ezxml_t superStruct, int subpool, const char *parentname, const char * corename);
int determine_struct_depth(int curLevel, ezxml_t superStruct);
int generate_struct_links(FILE *fd, int curLevel, ezxml_t superStruct, ezxml_t registry);
int generate_field_exchanges(FILE *fd, int curLevel, ezxml_t superStruct);
int generate_field_halo_exchanges_and_copies(ezxml_t registry);
int generate_field_inputs(FILE *fd, int curLevel, ezxml_t superStruct);
int generate_field_outputs(FILE *fd, int curLevel, ezxml_t superStruct);
int generate_field_reads_and_writes(ezxml_t registry);
int push_attributes(ezxml_t currentPosition);
int merge_structs_and_var_arrays(ezxml_t currentPosition);
int merge_streams(ezxml_t registry);
int parse_structs_from_registry(ezxml_t registry);

