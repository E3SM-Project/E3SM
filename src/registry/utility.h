// Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
// and the University Corporation for Atmospheric Research (UCAR).
//
// Unless noted otherwise source code is licensed under the BSD license.
// Additional copyright and license information can be found in the LICENSE file
// distributed with this code, or at http://mpas-dev.github.com/license.html
//


int is_derived_dim(char * d);
char * new_dimension_name(char * old_name);
void split_derived_dim_string(char * dim, char ** p1, char ** p2);
int is_integer_constant(char * c);
char * check_packages(ezxml_t registry, char * packages);
char * check_dimensions(ezxml_t registry, char * dims);
char * check_streams(char * streams);
int check_persistence(const char * persistence);
