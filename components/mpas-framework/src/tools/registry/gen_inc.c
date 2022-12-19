// Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
// and the University Corporation for Atmospheric Research (UCAR).
//
// Unless noted otherwise source code is licensed under the BSD license.
// Additional copyright and license information can be found in the LICENSE file
// distributed with this code, or at http://mpas-dev.github.io/license.html
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ezxml.h"
#include "registry_types.h"
#include "gen_inc.h"
#include "fortprintf.h"
#include "utility.h"

#define STR(s) #s
#define MACRO_TO_STR(s) STR(s)

void write_model_variables(ezxml_t registry){/*{{{*/
	const char * suffix = MACRO_TO_STR(MPAS_NAMELIST_SUFFIX);
	const char * exe_name = MACRO_TO_STR(MPAS_EXE_NAME);
	const char * git_ver = MACRO_TO_STR(MPAS_GIT_VERSION);

	const char *modelname, *corename, *version;
	FILE *fd;

	modelname = ezxml_attr(registry, "model");
	corename = ezxml_attr(registry, "core");
	version = ezxml_attr(registry, "version");

	fd = fopen("core_variables.inc", "w+");

	fortprintf(fd, "       core %% modelName = '%s'\n", modelname);
	fortprintf(fd, "       core %% coreName = '%s'\n", corename);
	fortprintf(fd, "       core %% modelVersion = '%s'\n", version);
	fortprintf(fd, "       core %% executableName = '%s'\n", exe_name);
	fortprintf(fd, "       core %% git_version = '%s'\n", git_ver);

	fclose(fd);

	fd = fopen("domain_variables.inc", "w+");

	fortprintf(fd, "       domain %% namelist_filename = 'namelist.%s'\n", suffix);
	fortprintf(fd, "       domain %% streams_filename = 'streams.%s'\n", suffix);

	fclose(fd);

}/*}}}*/


int write_field_pointer_arrays(FILE* fd){/*{{{*/
	fortprintf(fd, "\n");
	fortprintf(fd, "      type (field0DReal), pointer :: r0Ptr\n");
	fortprintf(fd, "      type (field1DReal), pointer :: r1Ptr\n");
	fortprintf(fd, "      type (field2DReal), pointer :: r2Ptr\n");
	fortprintf(fd, "      type (field3DReal), pointer :: r3Ptr\n");
	fortprintf(fd, "      type (field4DReal), pointer :: r4Ptr\n");
	fortprintf(fd, "      type (field5DReal), pointer :: r5Ptr\n");
	fortprintf(fd, "      type (field0DInteger), pointer :: i0Ptr\n");
	fortprintf(fd, "      type (field1DInteger), pointer :: i1Ptr\n");
	fortprintf(fd, "      type (field2DInteger), pointer :: i2Ptr\n");
	fortprintf(fd, "      type (field3DInteger), pointer :: i3Ptr\n");
	fortprintf(fd, "      type (field0DChar), pointer :: c0Ptr\n");
	fortprintf(fd, "      type (field1DChar), pointer :: c1Ptr\n");
	fortprintf(fd, "      type (field0DReal), dimension(:), pointer :: r0aPtr\n");
	fortprintf(fd, "      type (field1DReal), dimension(:), pointer :: r1aPtr\n");
	fortprintf(fd, "      type (field2DReal), dimension(:), pointer :: r2aPtr\n");
	fortprintf(fd, "      type (field3DReal), dimension(:), pointer :: r3aPtr\n");
	fortprintf(fd, "      type (field4DReal), dimension(:), pointer :: r4aPtr\n");
	fortprintf(fd, "      type (field5DReal), dimension(:), pointer :: r5aPtr\n");
	fortprintf(fd, "      type (field0DInteger), dimension(:), pointer :: i0aPtr\n");
	fortprintf(fd, "      type (field1DInteger), dimension(:), pointer :: i1aPtr\n");
	fortprintf(fd, "      type (field2DInteger), dimension(:), pointer :: i2aPtr\n");
	fortprintf(fd, "      type (field3DInteger), dimension(:), pointer :: i3aPtr\n");
	fortprintf(fd, "      type (field0DChar), dimension(:), pointer :: c0aPtr\n");
	fortprintf(fd, "      type (field1DChar), dimension(:), pointer :: c1aPtr\n");
	fortprintf(fd, "\n");

	return 0;
}/*}}}*/


int set_pointer_name(int type, int ndims, char *pointer_name, int time_levs){/*{{{*/

	char suffix[6];

	if (time_levs > 1) {
		snprintf(suffix, 6, "aPtr");
	} else {
		snprintf(suffix, 6, "Ptr");
	}

	if(type == REAL) {
		switch (ndims){
			default:
			case 0:
				snprintf(pointer_name, 1024, "r0%s", suffix);
				break;
			case 1:
				snprintf(pointer_name, 1024, "r1%s", suffix);
				break;
			case 2:
				snprintf(pointer_name, 1024, "r2%s", suffix);
				break;
			case 3:
				snprintf(pointer_name, 1024, "r3%s", suffix);
				break;
			case 4:
				snprintf(pointer_name, 1024, "r4%s", suffix);
				break;
			case 5:
				snprintf(pointer_name, 1024, "r5%s", suffix);
				break;
		}
	} else if (type == INTEGER) {
		switch (ndims){
			default:
			case 0:
				snprintf(pointer_name, 1024, "i0%s", suffix);
				break;
			case 1:
				snprintf(pointer_name, 1024, "i1%s", suffix);
				break;
			case 2:
				snprintf(pointer_name, 1024, "i2%s", suffix);
				break;
			case 3:
				snprintf(pointer_name, 1024, "i3%s", suffix);
				break;
		}
	} else if (type == CHARACTER) {
		switch (ndims){
			default:
			case 0:
				snprintf(pointer_name, 1024, "c0%s", suffix);
				break;
			case 1:
				snprintf(pointer_name, 1024, "c1%s", suffix);
				break;
		}
	}

	return 0;
}/*}}}*/


int add_package_to_list(const char * package, const char * package_list){/*{{{*/
	char *token, *string, *tofree;

	string = strdup(package_list);
	tofree = string;
	token = strsep(&string, ";");

	if(strcmp(package, token) == 0){
		return 0;
	}

	while( (token = strsep(&string, ";")) != NULL){
		if(strcmp(package, token) == 0){

			return 0;
		}
	}

	return 1;
}/*}}}*/


int build_struct_package_lists(ezxml_t currentPosition, char * out_packages){/*{{{*/
	ezxml_t child_xml1, child_xml2;

	const char *package_list;
	const char *name;

	char *token, *string, *tofree;
	int empty_packages;
	int empty_struct;

	package_list = ezxml_attr(currentPosition, "packages");

	empty_packages = 0;
	empty_struct = 1;

	// Check for vars that don't have packages.
	for(child_xml1 = ezxml_child(currentPosition, "var"); child_xml1 && !empty_packages; child_xml1 = child_xml1->next){
		package_list = ezxml_attr(child_xml1, "packages");

		if(!package_list){
			empty_packages = 1;
		}
		empty_struct = 0;
	}

	// Check for vararrays and constituents that don't have packages.
	for(child_xml1 = ezxml_child(currentPosition, "var_array"); child_xml1 && !empty_packages; child_xml1 = child_xml1->next){
		package_list = ezxml_attr(child_xml1, "packages");

		if(!package_list){
			for(child_xml2 = ezxml_child(child_xml1, "var"); child_xml2 && !empty_packages; child_xml2 = child_xml2->next){
				package_list = ezxml_attr(child_xml2, "packages");

				if(!package_list){
					empty_packages = 1;
				}
				empty_struct = 0;
			}
		}
	}

	// If any var/var_array doesn't have packages on it, the struct doesn't have packages on it.
	if(empty_packages || empty_struct){
		return 1;
	} else {
		// Build unique list of packages from nested vars and var arrays.
		for(child_xml1 = ezxml_child(currentPosition, "var_array"); child_xml1; child_xml1 = child_xml1->next){
			package_list = ezxml_attr(child_xml1, "packages");

			// Build list of unique packages from var_array
			if(package_list){
				string = strdup(package_list);
				tofree = string;
				token = strsep(&string, ";");

				if(out_packages[0] == '\0'){
					sprintf(out_packages, "%s", token);
				} else if(add_package_to_list(token, out_packages)){
					sprintf(out_packages, "%s;%s", out_packages, token);
				}

				while( (token = strsep(&string, ";")) != NULL){
					if(add_package_to_list(token, out_packages)){
						sprintf(out_packages, "%s;%s", out_packages, token);
					}
				}

				free(tofree);
			}

			for(child_xml2 = ezxml_child(child_xml1, "var"); child_xml2; child_xml2 = child_xml2->next){
				package_list = ezxml_attr(child_xml2, "packages");

				// Build list of unique packages from child var
				if(package_list){
					string = strdup(package_list);
					tofree = string;
					token = strsep(&string, ";");

					if(out_packages[0] == '\0'){
						sprintf(out_packages, "%s", token);
					} else if(add_package_to_list(token, out_packages)){
						sprintf(out_packages, "%s;%s", out_packages, token);
					}

					while( (token = strsep(&string, ";")) != NULL){
						if(add_package_to_list(token, out_packages)){
							sprintf(out_packages, "%s;%s", out_packages, token);
						}
					}

					free(tofree);
				}
			}
		}

		for(child_xml1 = ezxml_child(currentPosition, "var"); child_xml1; child_xml1 = child_xml1->next){
			package_list = ezxml_attr(child_xml1, "packages");

			// Build list of unique packages from child var
			if(package_list){
				string = strdup(package_list);
				tofree = string;
				token = strsep(&string, ";");

				if(out_packages[0] == '\0'){
					sprintf(out_packages, "%s", token);
				} else if(add_package_to_list(token, out_packages)){
					sprintf(out_packages, "%s;%s", out_packages, token);
				}

				while( (token = strsep(&string, ";")) != NULL){
					if(add_package_to_list(token, out_packages)){
						sprintf(out_packages, "%s;%s", out_packages, token);
					}
				}

				free(tofree);
			}
		}
		return 0;
	}
}/*}}}*/


int get_dimension_information(ezxml_t registry, const char *test_dimname, int *time_dim, int *decomp){/*{{{*/
	ezxml_t dims_xml, dim_xml;
	const char *dimname, *dimdecomp;

	int found;

	found = 0;
	(*time_dim) = 0;
	(*decomp) = -1;

	if ( strcmp(test_dimname, "Time") == 0 ){
		(*time_dim) = 1;
		found = 1;
		return found;
	}

	for (dims_xml = ezxml_child(registry, "dims"); dims_xml && !found; dims_xml = dims_xml->next){
		for (dim_xml = ezxml_child(dims_xml, "dim"); dim_xml && !found; dim_xml = dim_xml->next){
			dimname = ezxml_attr(dim_xml, "name");
			dimdecomp = ezxml_attr(dim_xml, "decomposition");

			if ( dimname != NULL && strcmp(dimname, test_dimname) == 0 ) {
				found = 1;

				/* Determine decomposition */
				if ( strcmp(dimname, "nCells") == 0 ) {
					if ( (*decomp) == -1 ) {
						(*decomp) = CELLS;
					} else {
						fprintf(stderr, "ERROR: Multiple decomposition types...\n");
						return 1;
					}
				} else if ( strcmp(dimname, "nEdges") == 0 ) {
					if ( (*decomp) == -1 ) {
						(*decomp) = EDGES;
					} else {
						fprintf(stderr, "ERROR: Multiple decomposition types...\n");
						return 1;
					}
				} else if ( strcmp(dimname, "nVertices") == 0 ) {
					if ( (*decomp) == -1 ) {
						(*decomp) = VERTICES;
					} else {
						fprintf(stderr, "ERROR: Multiple decomposition types...\n");
						return 1;
					}
				} else {
					if ( dimdecomp != NULL ) {
						if ( strcmp(dimdecomp, "none") != 0 ) {
							if ( (*decomp) == -1 ) {
							(*decomp) = INTERNAL;
					} else {
						fprintf(stderr, "ERROR: Multiple decomposition types...\n");
						return 1;
					}
						}
					}
				}
			}
		}
	}

	return found;
}/*}}}*/


int build_dimension_information(ezxml_t registry, ezxml_t var, int *ndims, int *has_time, int *decomp){/*{{{*/
	const char *vardims;
	const char *varname;

	int dim_exists, is_time, dim_decomp;

	char *string, *tofree, *token;

	varname = ezxml_attr(var, "name");
	vardims = ezxml_attr(var, "dimensions");
	(*ndims) = 0;
	(*decomp) = -1;
	(*has_time) = 0;

	string = strdup(vardims);
	tofree = string;
	token = strsep(&string, " ");

	/* Handle first dimension in list */
	if ( strlen(token) > 0 ) {
		dim_exists = get_dimension_information(registry, token, &is_time, &dim_decomp);

		if ( dim_exists && !is_time) {
			(*ndims)++;

			if ( (*decomp) == -1 ) {
				(*decomp) = dim_decomp;
			} else if ( (*decomp) != -1 && dim_decomp != -1) {
				fprintf(stderr, "ERROR: Variable %s contains multiple decomposed dimensions in list: %s\n", varname, vardims);
				return 1;
			}
        } else if (is_time) {
            (*has_time) = 1;
		} else if (!dim_exists) {
			fprintf(stderr, "ERROR: Dimension %s on variable %s doesn't exist.\n", token, varname);
			return 1;
		}
	}

	/* Handle remaining dimensions in list */
	while( (token = strsep(&string, " ")) != NULL ){
		dim_exists = get_dimension_information(registry, token, &is_time, &dim_decomp);

		if ( dim_exists && !is_time ) {
			(*ndims)++;

			if ( (*decomp) == -1 ) {
				(*decomp) = dim_decomp;
			} else if ( (*decomp) != -1 && dim_decomp != -1) {
				fprintf(stderr, "ERROR: Variable %s contains multiple decomposed dimensions in list: %s\n", varname, vardims);
				return 1;
			}
        } else if (is_time) {
            (*has_time) = 1;
		} else if (!dim_exists) {
			fprintf(stderr, "ERROR: Dimension %s on variable %s doesn't exist.\n", token, varname);
			return 1;
		}
	}

	free(tofree);

	return 0;
}/*}}}*/


int get_field_information(const char *vartype, const char *varval, char *default_value, const char *varmissval, char *missing_value, int *type){/*{{{*/
	if (strcmp(vartype, "real") == 0){
		(*type) = REAL;
		if(!varval){
			snprintf(default_value, 1024, "0.0");
		} else {
			snprintf(default_value, 1024, "%s", varval);
		}
	} else if (strcmp(vartype, "integer") == 0){
		(*type) = INTEGER;
		if(!varval){
			snprintf(default_value, 1024, "0");
		} else {
			snprintf(default_value, 1024, "%s", varval);
		}
	} else if (strcmp(vartype, "text") == 0){
		(*type) = CHARACTER;
		if(!varval){
			snprintf(default_value, 1024, "''");
		} else {
			snprintf(default_value, 1024, "%s", varval);
		}
	}

	if (strcmp(vartype, "real") == 0){
		(*type) = REAL;
		if(!varmissval || strcmp(varmissval, "FILLVAL") == 0){
			snprintf(missing_value, 1024, "MPAS_REAL_FILLVAL");
		} else {
			snprintf(missing_value, 1024, "%s", varmissval);
		}
	} else if (strcmp(vartype, "integer") == 0){
		(*type) = INTEGER;
		if(!varmissval || strcmp(varmissval, "FILLVAL") == 0){
			snprintf(missing_value, 1024, "MPAS_INT_FILLVAL");
		} else {
			snprintf(missing_value, 1024, "%s", varmissval);
		}
	} else if (strcmp(vartype, "text") == 0){
		(*type) = CHARACTER;
		if(!varmissval || strcmp(varmissval, "FILLVAL") == 0){
			snprintf(missing_value, 1024, "MPAS_CHAR_FILLVAL");
		} else {
			snprintf(missing_value, 1024, "'%s'", varmissval);
		}
	}

	return 0;
}/*}}}*/


int parse_packages_from_registry(ezxml_t registry)/*{{{*/
{
	ezxml_t packages_xml, package_xml;

	const char *packagename, *packagedesc, *const_core;
	FILE *fd;
	char core_string[1024];

	const_core = ezxml_attr(registry, "core_abbrev");

	fd = fopen("define_packages.inc", "w+");

	sprintf(core_string, "%s", const_core);

	fortprintf(fd, "   function %s_define_packages(packagePool) result(iErr)\n", core_string);
	fortprintf(fd, "      use mpas_derived_types\n");
	fortprintf(fd, "      use mpas_pool_routines\n");
	fortprintf(fd, "      use mpas_io_units\n");
	fortprintf(fd, "      implicit none\n");
	fortprintf(fd, "      type (mpas_pool_type), intent(inout) :: packagePool !< Input: MPAS Pool for containing package logicals.\n\n");
	fortprintf(fd, "      integer :: iErr\n");
	fortprintf(fd, "\n");
	fortprintf(fd, "      iErr = 0\n");

	// Parse Packages
	for (packages_xml = ezxml_child(registry, "packages"); packages_xml; packages_xml = packages_xml->next){
		for (package_xml = ezxml_child(packages_xml, "package"); package_xml; package_xml = package_xml->next){
			packagename = ezxml_attr(package_xml, "name");
			packagedesc = ezxml_attr(package_xml, "description");

			fortprintf(fd, "      call mpas_pool_add_package(packagePool, '%sActive', .false.)\n", packagename);
		}
	}

	fortprintf(fd, "   end function %s_define_packages\n", core_string);

	fclose(fd);

	return 0;
}/*}}}*/


int parse_namelist_records_from_registry(ezxml_t registry)/*{{{*/
{
	ezxml_t nmlrecs_xml, nmlopt_xml;

	const char *const_core;
	const char *nmlrecname, *nmlrecindef, *nmlrecinsub;
	const char *nmloptname, *nmlopttype, *nmloptval, *nmloptunits, *nmloptdesc, *nmloptposvals, *nmloptindef;
	const char *nmloptbounds, *nmloptcellmeasures, *nmloptcellmethod, *nmloptcoords, *nmloptstdname;

	char pool_name[1024];
	char core_string[1024];

	int in_subpool;

	FILE *fd, *fd2, *fcd, *fcg;

	const_core = ezxml_attr(registry, "core_abbrev");

	sprintf(core_string, "%s", const_core);

	fd = fopen("namelist_defines.inc", "w+");
	fd2 = fopen("namelist_call.inc", "w+");
	fcd = fopen("config_declare.inc", "w+");
	fcg = fopen("config_get.inc", "w+");

	fortprintf(fd2, "   function %s_setup_namelists(configPool, namelistFilename, dminfo) result(iErr)\n", core_string);
	fortprintf(fd2, "      use mpas_derived_types\n");
	fortprintf(fd2, "      use mpas_pool_routines\n");
	fortprintf(fd2, "      use mpas_io_units\n");
	fortprintf(fd2, "      use mpas_abort, only : mpas_dmpar_global_abort\n");
	fortprintf(fd2, "      use mpas_log, only : mpas_log_write\n");
	fortprintf(fd2, "      implicit none\n");
	fortprintf(fd2, "      type (mpas_pool_type), intent(inout) :: configPool\n");
	fortprintf(fd2, "      character (len=*), intent(in) :: namelistFilename\n");
	fortprintf(fd2, "      type (dm_info), intent(in) :: dminfo\n");
	fortprintf(fd2, "      integer :: iErr\n");
	fortprintf(fd2, "\n");
	fortprintf(fd2, "      integer :: unitNumber\n");
	fortprintf(fd2, "      logical :: nmlExists\n");
	fortprintf(fd2, "\n");
	fortprintf(fd2, "      iErr = 0\n");
	fortprintf(fd2, "      unitNumber = 21\n");
	fortprintf(fd2, "      call mpas_log_write('Reading namelist from file '//trim(namelistFilename))\n");
	fortprintf(fd2, "      inquire(file=trim(namelistFilename), exist=nmlExists)\n");
	fortprintf(fd2, "      if ( .not. nmlExists ) then\n");
	fortprintf(fd2, "         call mpas_dmpar_global_abort('ERROR: Namelist file '//trim(namelistFilename)//' does not exist.')\n");
	fortprintf(fd2, "      end if\n");
	fortprintf(fd2, "      open(unitNumber,file=trim(namelistFilename),status='old',form='formatted')\n");
	fortprintf(fd2, "\n");


	// Parse Namelist Records
	for (nmlrecs_xml = ezxml_child(registry, "nml_record"); nmlrecs_xml; nmlrecs_xml = nmlrecs_xml->next){
		nmlrecname = ezxml_attr(nmlrecs_xml, "name");
		nmlrecindef = ezxml_attr(nmlrecs_xml, "in_defaults");
		nmlrecinsub = ezxml_attr(nmlrecs_xml, "in_subpool");

		in_subpool = 0;

		if(nmlrecinsub){
			if(strcmp(nmlrecinsub, "true") == 0){
				in_subpool = 1;
			}
		}

		if(in_subpool){
			sprintf(pool_name, "recordPool");
		} else {
			sprintf(pool_name, "configPool");
		}

		// Add call to driver routine.
		fortprintf(fd2, "      call %s_setup_nmlrec_%s(configPool, unitNumber, dminfo)\n", core_string, nmlrecname);

		// Start defining new subroutine for namelist record.
		fortprintf(fd, "   subroutine %s_setup_nmlrec_%s(configPool, unitNumber, dminfo)\n", core_string, nmlrecname);
		fortprintf(fd, "      use mpas_log, only : mpas_log_write, mpas_log_escape_dollars\n");
		fortprintf(fd, "      implicit none\n");
		fortprintf(fd, "      type (mpas_pool_type), intent(inout) :: configPool\n");
		fortprintf(fd, "      integer, intent(in) :: unitNumber\n");
		fortprintf(fd, "      type (dm_info), intent(in) :: dminfo\n");
		fortprintf(fd, "      type (mpas_pool_type), pointer :: recordPool\n");
		fortprintf(fd, "      integer :: ierr\n");
		fortprintf(fd, "\n");

		// Define variable definitions prior to reading the namelist in.
		for (nmlopt_xml = ezxml_child(nmlrecs_xml, "nml_option"); nmlopt_xml; nmlopt_xml = nmlopt_xml->next){
			nmloptname = ezxml_attr(nmlopt_xml, "name");
			nmlopttype = ezxml_attr(nmlopt_xml, "type");
			nmloptval = ezxml_attr(nmlopt_xml, "default_value");
			nmloptunits = ezxml_attr(nmlopt_xml, "units");
			nmloptbounds = ezxml_attr(nmlopt_xml, "bounds");
			nmloptcellmeasures = ezxml_attr(nmlopt_xml, "cell_measures");
			nmloptcellmethod = ezxml_attr(nmlopt_xml, "cell_method");
			nmloptcoords = ezxml_attr(nmlopt_xml, "coordinates");
			nmloptstdname = ezxml_attr(nmlopt_xml, "standard_name");
			nmloptdesc = ezxml_attr(nmlopt_xml, "description");
			nmloptposvals = ezxml_attr(nmlopt_xml, "possible_values");
			nmloptindef = ezxml_attr(nmlopt_xml, "in_defaults");

			if(strncmp(nmlopttype, "real", 1024) == 0){
				fortprintf(fd, "      real (kind=RKIND) :: %s = %lf\n", nmloptname, (double)atof(nmloptval));
				fortprintf(fcd, "      real (kind=RKIND) :: %s\n", nmloptname);
			} else if(strncmp(nmlopttype, "integer", 1024) == 0){
				fortprintf(fd, "      integer :: %s = %d\n", nmloptname, atoi(nmloptval));
				fortprintf(fcd, "      integer :: %s\n", nmloptname);
			} else if(strncmp(nmlopttype, "logical", 1024) == 0){
				fortprintf(fcd, "      logical :: %s\n", nmloptname);
				if(strncmp(nmloptval, "true", 1024) == 0 || strncmp(nmloptval, ".true.", 1024) == 0){
					fortprintf(fd, "      logical :: %s = .true.\n", nmloptname);
				} else {
					fortprintf(fd, "      logical :: %s = .false.\n", nmloptname);
				}
			} else if(strncmp(nmlopttype, "character", 1024) == 0){
					fortprintf(fd, "      character (len=StrKIND) :: %s = '%s'\n", nmloptname, nmloptval);
					fortprintf(fcd, "      character (len=StrKIND) :: %s\n", nmloptname);
			}
		}
		fortprintf(fd, "\n");
		fortprintf(fcd, "\n");

		// Define the namelist block, to read the namelist record in.
		fortprintf(fd, "      namelist /%s/ &\n", nmlrecname);
		for (nmlopt_xml = ezxml_child(nmlrecs_xml, "nml_option"); nmlopt_xml; nmlopt_xml = nmlopt_xml->next){
			nmloptname = ezxml_attr(nmlopt_xml, "name");
			if(nmlopt_xml->next){
				fortprintf(fd, "         %s, &\n", nmloptname);
			} else {
				fortprintf(fd, "         %s\n", nmloptname);
			}
		}

		if(in_subpool){
			fortprintf(fd, "\n");
			fortprintf(fd, "      allocate(recordPool)\n");
			fortprintf(fd, "      call mpas_pool_create_pool(recordPool)\n");
			fortprintf(fd, "      call mpas_pool_add_subpool(configPool, '%s', recordPool)\n", nmlrecname);
			fortprintf(fd, "\n");
		}

		fortprintf(fd, "      if (dminfo %% my_proc_id == IO_NODE) then\n");
		fortprintf(fd, "! Rewinding before each read leads to errors when the code is built with\n");
		fortprintf(fd, "! the NAG Fortran compiler. If building with NAG, be kind and don't rewind.\n");
		fortprintf(fd, "#ifndef NAG_COMPILER\n");
		fortprintf(fd, "         rewind(unitNumber)\n");
		fortprintf(fd, "#endif\n");
		fortprintf(fd, "         read(unitNumber, %s, iostat=ierr)\n", nmlrecname);
		fortprintf(fd, "      end if\n");

		// Broadcast ierr, to check if a broadcast should happen for the options (if namelist was read in)
		fortprintf(fd, "      call mpas_dmpar_bcast_int(dminfo, ierr)\n");

		fortprintf(fd, "\n");
		// Define broadcast calls for namelist values.
		fortprintf(fd, "      if (ierr <= 0) then\n");
		for (nmlopt_xml = ezxml_child(nmlrecs_xml, "nml_option"); nmlopt_xml; nmlopt_xml = nmlopt_xml->next){
			nmloptname = ezxml_attr(nmlopt_xml, "name");
			nmlopttype = ezxml_attr(nmlopt_xml, "type");

			if(strncmp(nmlopttype, "real", 1024) == 0){
				fortprintf(fd, "         call mpas_dmpar_bcast_real(dminfo, %s)\n", nmloptname);
			} else if(strncmp(nmlopttype, "integer", 1024) == 0){
				fortprintf(fd, "         call mpas_dmpar_bcast_int(dminfo, %s)\n", nmloptname);
			} else if(strncmp(nmlopttype, "logical", 1024) == 0){
				fortprintf(fd, "         call mpas_dmpar_bcast_logical(dminfo, %s)\n", nmloptname);
			} else if(strncmp(nmlopttype, "character", 1024) == 0){
				fortprintf(fd, "         call mpas_dmpar_bcast_char(dminfo, %s)\n", nmloptname);
			}
		}
		fortprintf(fd, "         if (ierr < 0) then\n");
		fortprintf(fd, "            call mpas_log_write('*** Encountered an issue while attempting to read namelist record %s')\n", nmlrecname);
		fortprintf(fd, "            call mpas_log_write('    The following values will be used for variables in this record:')\n");
		fortprintf(fd, "            call mpas_log_write(' ')\n");
		for (nmlopt_xml = ezxml_child(nmlrecs_xml, "nml_option"); nmlopt_xml; nmlopt_xml = nmlopt_xml->next){
			nmloptname = ezxml_attr(nmlopt_xml, "name");
			nmlopttype = ezxml_attr(nmlopt_xml, "type");

			if (strncmp(nmlopttype, "character", 1024) == 0) {
				fortprintf(fd, "            call mpas_log_write('        %s = '//mpas_log_escape_dollars(%s))\n", nmloptname, nmloptname);
			}
			else if (strncmp(nmlopttype, "integer", 1024) == 0) {
				fortprintf(fd, "            call mpas_log_write('        %s = $i', intArgs=(/%s/))\n", nmloptname, nmloptname);
			}
			else if (strncmp(nmlopttype, "real", 1024) == 0) {
				fortprintf(fd, "            call mpas_log_write('        %s = $r', realArgs=(/%s/))\n", nmloptname, nmloptname);
			}
			else if (strncmp(nmlopttype, "logical", 1024) == 0) {
				fortprintf(fd, "            call mpas_log_write('        %s = $l', logicArgs=(/%s/))\n", nmloptname, nmloptname);
			}
			else {
				fortprintf(fd, "            call mpas_log_write('        %s = <unhandled namelist type>')\n", nmloptname);
			}
		}
		fortprintf(fd, "            call mpas_log_write(' ')\n");
		fortprintf(fd, "         end if\n");
		fortprintf(fd, "      else if (ierr > 0) then\n");
		fortprintf(fd, "         call mpas_log_write('Error while reading namelist record %s.', MPAS_LOG_CRIT)\n", nmlrecname);
		fortprintf(fd, "      end if\n");
		fortprintf(fd, "\n");

		for (nmlopt_xml = ezxml_child(nmlrecs_xml, "nml_option"); nmlopt_xml; nmlopt_xml = nmlopt_xml->next){
			nmloptname = ezxml_attr(nmlopt_xml, "name");

			fortprintf(fd, "      call mpas_pool_add_config(%s, '%s', %s)\n", pool_name, nmloptname, nmloptname);
			fortprintf(fcg, "      call mpas_pool_get_config_scalar(configPool, '%s', %s)\n", nmloptname, nmloptname);
		}
		fortprintf(fd, "\n");
		fortprintf(fcg, "\n");

		// End new subroutine for namelist record.
		fortprintf(fd, "   end subroutine %s_setup_nmlrec_%s\n", core_string, nmlrecname);
		fortprintf(fd, "\n\n");
	}

	fortprintf(fd2, "\n");
	fortprintf(fd2, "      close(unitNumber)\n");
	fortprintf(fd2, "   end function %s_setup_namelists\n", core_string);

	fclose(fd);
	fclose(fd2);
	fclose(fcd);
	fclose(fcg);

	return 0;
}/*}}}*/


int parse_dimensions_from_registry(ezxml_t registry)/*{{{*/
{
	ezxml_t dims_xml, dim_xml;
	ezxml_t nmlrec_xml, nmlopt_xml;

	const char *nmlrecname, *nmlrecinsub, *nmloptname, *nmlopttype;
	const char *dimname, *dimunits, *dimdesc, *dimdef, *dimdecomp;
	const char *corename;

	char option_name[1024];
	char core_string[1024];

	FILE *fd;

	int in_subpool, readable;

	corename = ezxml_attr(registry, "core_abbrev");

	sprintf(core_string, "%s", corename);

	fd = fopen("block_dimension_routines.inc", "w+");

	fortprintf(fd, "   function %s_setup_derived_dimensions(readDimensions, dimensionPool, configPool) result(iErr)\n", core_string);
	fortprintf(fd, "\n");
	fortprintf(fd, "      use mpas_derived_types\n");
	fortprintf(fd, "      use mpas_pool_routines\n");
	fortprintf(fd, "      use mpas_io_units\n");
	fortprintf(fd, "      use mpas_log, only : mpas_log_write\n");
	fortprintf(fd, "\n");
	fortprintf(fd, "      implicit none\n");
	fortprintf(fd, "\n");
	fortprintf(fd, "      type (mpas_pool_type), intent(inout) :: readDimensions !< Input: Pool to pull read dimensions from\n");
	fortprintf(fd, "      type (mpas_pool_type), intent(inout) :: configPool !< Input: Pool containing namelist options with configs\n");
	fortprintf(fd, "      type (mpas_pool_type), intent(inout) :: dimensionPool !< Input/Output: Pool to add dimensions into\n");
	fortprintf(fd, "\n");
	fortprintf(fd, "      integer :: iErr, errLevel\n");
	fortprintf(fd, "\n");

	// Define all dimensions
	for (dims_xml = ezxml_child(registry, "dims"); dims_xml; dims_xml = dims_xml->next){
		for (dim_xml = ezxml_child(dims_xml, "dim"); dim_xml; dim_xml = dim_xml->next){
			dimname = ezxml_attr(dim_xml, "name");

			fortprintf(fd, "      integer, pointer :: %s\n", dimname);
		}
	}


	// Define all namelist options that will be used
	for (dims_xml = ezxml_child(registry, "dims"); dims_xml; dims_xml = dims_xml->next){
		for (dim_xml = ezxml_child(dims_xml, "dim"); dim_xml; dim_xml = dim_xml->next){
			dimname = ezxml_attr(dim_xml, "name");
			dimdef = ezxml_attr(dim_xml, "definition");

			if(dimdef != NULL){
				/* Namelist defined dimension */
				if(strncmp(dimdef, "namelist:", 9) == 0){
					snprintf(option_name, 1024, "%s", (dimdef)+9);
					/* Need to define a variable to hold the namelist value */
					/* First need to find the registry defined namlist option, so we can determine type: */
					for (nmlrec_xml = ezxml_child(registry, "nml_record"); nmlrec_xml; nmlrec_xml = nmlrec_xml->next){
						nmlrecname = ezxml_attr(nmlrec_xml, "name");
						for (nmlopt_xml = ezxml_child(nmlrec_xml, "nml_option"); nmlopt_xml; nmlopt_xml = nmlopt_xml->next){
							nmloptname = ezxml_attr(nmlopt_xml, "name");
							nmlopttype = ezxml_attr(nmlopt_xml, "type");

							if(strncmp(option_name, nmloptname, 1024) == 0){
								if(strncmp(nmlopttype, "real", 1024) == 0){
									fortprintf(fd, "      real (kind=RKIND), pointer :: %s\n", nmloptname);
								} else if(strncmp(nmlopttype, "integer", 1024) == 0){
									fortprintf(fd, "      integer, pointer :: %s\n", nmloptname);
								} else if(strncmp(nmlopttype, "logical", 1024) == 0){
									fortprintf(fd, "      logical, pointer :: %s\n", nmloptname);
								} else if(strncmp(nmlopttype, "character", 1024) == 0){
									fortprintf(fd, "      character (len=StrKIND), pointer :: %s\n", nmloptname);
								}
							}
						}
					}
				}
			}
		}
	}

	fortprintf(fd, "\n");

	fortprintf(fd, "      iErr = 0\n");
	fortprintf(fd, "      errLevel = mpas_pool_get_error_level()\n");
	fortprintf(fd, "      call mpas_pool_set_error_level(MPAS_POOL_SILENT)\n");

	fortprintf(fd, "\n");

	/* Get all namelist options first, so any derived dimensions can be properly defined based on them */
	for (dims_xml = ezxml_child(registry, "dims"); dims_xml; dims_xml = dims_xml->next){
		for (dim_xml = ezxml_child(dims_xml, "dim"); dim_xml; dim_xml = dim_xml->next){
			dimname = ezxml_attr(dim_xml, "name");
			dimdef = ezxml_attr(dim_xml, "definition");

			if(dimdef != NULL){
				/* Namelist defined dimension */
				if(strncmp(dimdef, "namelist:", 9) == 0){
					snprintf(option_name, 1024, "%s", (dimdef)+9);

					fortprintf(fd, "      nullify(%s)\n", option_name);
					fortprintf(fd, "      call mpas_pool_get_config(configPool, '%s', %s)\n", option_name, option_name);
				}
			}
		}
	}

	fortprintf(fd, "\n");

	/* Get all dimensions first, so any derived dimensions can be properly defined based on them */
	for (dims_xml = ezxml_child(registry, "dims"); dims_xml; dims_xml = dims_xml->next) {
		for (dim_xml = ezxml_child(dims_xml, "dim"); dim_xml; dim_xml = dim_xml->next) {
			dimname = ezxml_attr(dim_xml, "name");

			fortprintf(fd, "      nullify(%s)\n", dimname);
			fortprintf(fd, "      call mpas_pool_get_dimension(dimensionPool, '%s', %s)\n", dimname, dimname);
		}
	}

	fortprintf(fd, "\n");

	fortprintf(fd, "call mpas_log_write('Assigning remaining dimensions from definitions in Registry.xml ...')\n");

	for (dims_xml = ezxml_child(registry, "dims"); dims_xml; dims_xml = dims_xml->next) {
		for (dim_xml = ezxml_child(dims_xml, "dim"); dim_xml; dim_xml = dim_xml->next) {
			dimname = ezxml_attr(dim_xml, "name");
			dimdef = ezxml_attr(dim_xml, "definition");

			/* If dimension has a definition, check if the value of the dim is NaN, then write the definition */
			if ( dimdef != NULL ) {
				fortprintf(fd, "      call mpas_pool_get_dimension(dimensionPool, '%s', %s)\n", dimname, dimname);
				fortprintf(fd, "      if ( .not. associated(%s) ) then\n", dimname);
				fortprintf(fd, "         allocate(%s)\n", dimname);
				// Namelist defined dimension
				if(strncmp(dimdef, "namelist:", 9) == 0){
					snprintf(option_name, 1024, "%s", (dimdef)+9);
					fortprintf(fd, "         %s = %s\n", dimname, option_name);
					fortprintf(fd, "call mpas_log_write('       %s = $i (%s)', intArgs=(/%s/))\n", dimname, option_name, option_name);
				} else {
					fortprintf(fd, "         %s = %s\n", dimname, dimdef);
					fortprintf(fd, "call mpas_log_write('       %s = $i', intArgs=(/%s/))\n", dimname, dimdef);
				}
				fortprintf(fd, "         call mpas_pool_add_dimension(dimensionPool, '%s', %s)\n", dimname, dimname);

				fortprintf(fd, "          else if ( %s == MPAS_MISSING_DIM ) then\n", dimname, dimname);
				// Namelist defined dimension
				if(strncmp(dimdef, "namelist:", 9) == 0){
					snprintf(option_name, 1024, "%s", (dimdef)+9);
					fortprintf(fd, "         %s = %s\n", dimname, option_name);
				} else {
					fortprintf(fd, "         %s = %s\n", dimname, dimdef);
				}

				fortprintf(fd, "          end if\n\n");
			} else {
				fortprintf(fd, "      if ( .not. associated(%s) ) then\n", dimname);
				fortprintf(fd, "         allocate(%s)\n", dimname);
				fortprintf(fd, "         %s = MPAS_MISSING_DIM\n", dimname);
				fortprintf(fd, "         call mpas_pool_add_dimension(dimensionPool, '%s', %s)\n", dimname, dimname);
				fortprintf(fd, "      end if\n\n");
			}
		}
	}
	fortprintf(fd, "      call mpas_log_write(' ')\n");
	fortprintf(fd, "      call mpas_log_write(' ----- done assigning dimensions from Registry.xml -----')\n");
	fortprintf(fd, "      call mpas_log_write(' ')\n");
	fortprintf(fd, "      call mpas_log_write(' ')\n");

	fortprintf(fd, "      call mpas_pool_set_error_level(errLevel)\n\n");

	fortprintf(fd, "   end function %s_setup_derived_dimensions\n", core_string);

	fortprintf(fd, "\n\n");

	fortprintf(fd, "   function %s_setup_decomposed_dimensions(block, manager, readDimensions, dimensionPool, totalBlocks) result(iErr)\n", core_string);
	fortprintf(fd, "\n");
	fortprintf(fd, "      use mpas_derived_types\n");
	fortprintf(fd, "      use mpas_decomp\n");
	fortprintf(fd, "      use mpas_pool_routines\n");
	fortprintf(fd, "      use mpas_io_units\n");
	fortprintf(fd, "      use mpas_abort, only : mpas_dmpar_global_abort\n");
	fortprintf(fd, "      use mpas_log, only : mpas_log_write\n");
	fortprintf(fd, "\n");
	fortprintf(fd, "      implicit none\n");
	fortprintf(fd, "\n");
	fortprintf(fd, "      type (block_type), intent(inout) :: block !< Input: Pointer to block\n");
	fortprintf(fd, "      type (mpas_streamManager_type), intent(inout) :: manager !< Input: Stream manager\n");
	fortprintf(fd, "      type (mpas_pool_type), intent(inout) :: readDimensions !< Input: Pool to pull read dimensions from\n");
	fortprintf(fd, "      type (mpas_pool_type), intent(inout) :: dimensionPool !< Input/Output: Pool to add dimensions into\n");
	fortprintf(fd, "      integer, intent(in) :: totalBlocks !< Input: Number of blocks\n");
	fortprintf(fd, "\n");
	fortprintf(fd, "      integer :: iErr\n");
	fortprintf(fd, "      type (field1DInteger), pointer :: ownedIndices\n");
	fortprintf(fd, "      procedure (mpas_decomp_function), pointer :: decompFunc\n");
	fortprintf(fd, "\n");

	/* Define decomposed dimension integers */
	for (dims_xml = ezxml_child(registry, "dims"); dims_xml; dims_xml = dims_xml->next) {
		for (dim_xml = ezxml_child(dims_xml, "dim"); dim_xml; dim_xml = dim_xml->next) {
			dimname = ezxml_attr(dim_xml, "name");
			dimdecomp = ezxml_attr(dim_xml, "decomposition");

			if ( dimdecomp != NULL && strcmp(dimdecomp, "none") != 0 ) {
				fortprintf(fd, "      integer, pointer :: %s\n", dimname);
			}
		}
	}

	fortprintf(fd, "\n");
	fortprintf(fd, "      iErr = 0\n");
	fortprintf(fd, "      call mpas_log_write('Processing decomposed dimensions ...')\n\n");

	/* Retrieve dimension integers */
	for (dims_xml = ezxml_child(registry, "dims"); dims_xml; dims_xml = dims_xml->next) {
		for (dim_xml = ezxml_child(dims_xml, "dim"); dim_xml; dim_xml = dim_xml->next) {
			dimname = ezxml_attr(dim_xml, "name");
			dimdecomp = ezxml_attr(dim_xml, "decomposition");

			if ( dimdecomp != NULL && strcmp(dimdecomp, "none") != 0 ) {
				fortprintf(fd, "      call mpas_pool_get_dimension(readDimensions, '%s', %s)\n", dimname, dimname);
				fortprintf(fd, "      if ( .not. associated(%s)) then\n", dimname);
				fortprintf(fd, "         call mpas_log_write('Dimension ''%s'' was not defined, and cannot be decomposed.', MPAS_LOG_WARN)\n", dimname);
				fortprintf(fd, "      else\n");
				fortprintf(fd, "         call mpas_decomp_get_method(block %% domain %% decompositions, '%s', decompFunc, iErr)\n", dimdecomp);
				fortprintf(fd, "         if ( iErr /= MPAS_DECOMP_NOERR ) then\n");
				fortprintf(fd, "            call mpas_dmpar_global_abort('ERROR: Decomposition method \'\'%s\'\' used by dimension \'\'%s\'\' does not exist.')\n", dimdecomp, dimname);
				fortprintf(fd, "         end if\n");
				fortprintf(fd, "\n");
				fortprintf(fd, "         allocate(ownedIndices)\n");
				fortprintf(fd, "         ownedIndices %% hasTimeDimension = .false.\n");
				fortprintf(fd, "         ownedIndices %% isActive = .true.\n");
				fortprintf(fd, "         ownedIndices %% isVarArray = .false.\n");
				fortprintf(fd, "         ownedIndices %% isDecomposed = .false.\n");
				fortprintf(fd, "         ownedIndices %% isPersistent = .true.\n");
				fortprintf(fd, "         ownedIndices %% defaultValue = 0\n");
				fortprintf(fd, "         ownedIndices %% fieldName = '%sOwnedIndices'\n", dimname);
				fortprintf(fd, "         ownedIndices %% dimNames(1) = '%s'\n", dimname);
				fortprintf(fd, "         iErr = decompFunc(block, manager, %s, totalBlocks, ownedIndices %% array)\n", dimname);
				fortprintf(fd, "         ownedIndices %% dimSizes(1) = size(ownedIndices %% array, dim=1)\n");
				fortprintf(fd, "         call mpas_pool_add_field(block %% allFields, '%sOwnedIndices', ownedIndices)\n", dimname);
				fortprintf(fd, "         call mpas_pool_get_dimension(block %% dimensions, '%s', %s)\n", dimname, dimname);
				fortprintf(fd, "         %s = size(ownedIndices %% array, dim=1)\n", dimname);
				fortprintf(fd, "         call mpas_log_write('       %s => $i indices owned by block $i', intArgs=(/%s, block %% blockID/))\n", dimname, dimname);
				fortprintf(fd, "      end if\n");
				fortprintf(fd, "\n");
			}
		}
	}

	fortprintf(fd, "      call mpas_log_write(' ')\n");
	fortprintf(fd, "      call mpas_log_write(' ----- done processing decomposed dimensions -----')\n");
	fortprintf(fd, "      call mpas_log_write(' ')\n");
	fortprintf(fd, "      call mpas_log_write(' ')\n");

	fortprintf(fd, "\n");
	fortprintf(fd, "   end function %s_setup_decomposed_dimensions\n", core_string);

	fclose(fd);

	return 0;
}/*}}}*/


int parse_var_array(FILE *fd, ezxml_t registry, ezxml_t superStruct, ezxml_t varArray, const char * corename)/*{{{*/
{
	ezxml_t struct_xml, var_arr_xml, var_xml, var_xml2;
	ezxml_t packages_xml, package_xml;
	ezxml_t streams_xml, stream_xml, streams_xml2, stream_xml2;

	const char *structname, *structlevs, *structpackages;
	const char *substructname;
	const char *vararrname, *vararrtype, *vararrdims, *vararrpersistence, *vararrdefaultval, *vararrpackages, *vararrmissingval, *vararrmissing_value_mask;
	const char *vararrbounds, *vararrcellmeasures, *vararrcellmethod, *vararrcoords, *vararrstdname;
	const char *varname, *varpersistence, *vartype, *vardims, *varunits, *vardesc, *vararrgroup, *varstreams, *vardefaultval, *varpackages;
	const char *varbounds, *varcellmeasures, *varcellmethod, *varcoords, *varstdname;
	const char *varname2, *vararrgroup2, *vararrname_in_code;
	const char *varname_in_code;
	const char *varname_in_output;
	const char *streamname, *streamname2;
	const char *packagename;
	const char *vararrtimelevs;

	char package_list[2048];
	int no_packages;

	int err;

	int iostreams;
	int time_lev, time_levs;
	int i, skip_var, skip_stream;
	int ndims, type, hasTime, decomp, in_stream;
	int persistence;
	char *string, *tofree, *token;
	char temp_str[1024];
	char pointer_name[1024];
	char pointer_name_arr[1024];
	char spacing[1024], sub_spacing[1024];
	char default_value[1024];
	char missing_value[1024];

	structname = ezxml_attr(superStruct, "name");

	var_arr_xml = varArray;

	// All sub-structs have been parsed and generated at this point. Time to generate this struct
	// Start by generating variable arrays
	vararrname = ezxml_attr(var_arr_xml, "name");
	vararrtype = ezxml_attr(var_arr_xml, "type");
	vararrdims = ezxml_attr(var_arr_xml, "dimensions");
	vararrpersistence = ezxml_attr(var_arr_xml, "persistence");
	vararrdefaultval = ezxml_attr(var_arr_xml, "default_value");
	vararrmissingval = ezxml_attr(var_arr_xml, "missing_value");
	vararrpackages = ezxml_attr(var_arr_xml, "packages");
	vararrtimelevs = ezxml_attr(var_arr_xml, "time_levs");
	vararrname_in_code = ezxml_attr(var_arr_xml, "name_in_code");
	vararrmissing_value_mask= ezxml_attr(var_arr_xml, "missing_value_mask");
	vararrbounds = ezxml_attr(var_arr_xml, "bounds");
	vararrcellmeasures = ezxml_attr(var_arr_xml, "cell_measures");
	vararrcellmethod = ezxml_attr(var_arr_xml, "cell_method");
	vararrcoords = ezxml_attr(var_arr_xml, "coordinates");
	vararrstdname = ezxml_attr(var_arr_xml, "standard_name");

	package_list[0] = '\0';
	no_packages = build_struct_package_lists(varArray, package_list);

	if(!vararrtimelevs){
		vararrtimelevs = ezxml_attr(superStruct, "time_levs");
	}

	if(vararrtimelevs){
		time_levs = atoi(vararrtimelevs);
		if(time_levs < 1){
			time_levs = 1;
		}
	} else {
		time_levs = 1;
	}

	if(!vararrname_in_code){
		vararrname_in_code = ezxml_attr(var_arr_xml, "name");
	}

	persistence = check_persistence(vararrpersistence);

	fortprintf(fd, "! Define var array %s\n", vararrname);

	snprintf(spacing, 1024, "      ");

	// Determine field type and default value.
	get_field_information(vararrtype, vararrdefaultval, default_value, vararrmissingval, missing_value, &type);

	// If a default_value is not specified, but a missing_value is, then set the
	// default_value (defaultValue) to missing_value
	if(!vararrdefaultval && vararrmissingval) {
		snprintf(default_value, 1024, "%s ! defaultValue taking specified missing_value", missing_value);
	}

	// Determine ndims, hasTime, and decomp type
	build_dimension_information(registry, var_arr_xml, &ndims, &hasTime, &decomp);
	ndims++; // Add a dimension for constituents in var_array

	// Determine name of pointer for this field.
	set_pointer_name(type, ndims, pointer_name, time_levs);
	if (time_levs > 1) {
		fortprintf(fd, "      allocate(%s(%d))\n", pointer_name, time_levs);
	} else {
		fortprintf(fd, "      allocate(%s)\n", pointer_name);
	}

	fortprintf(fd, "      index_counter = 0\n", spacing);
	fortprintf(fd, "      group_counter = -1\n", spacing);
	fortprintf(fd, "      group_start = -1\n", spacing);
	fortprintf(fd, "      group_started = .false.\n", spacing);
	fortprintf(fd, "\n");

	// Write index values and group counter values.
	// Define each array_group in contiguous sections.
	for (var_xml = ezxml_child(var_arr_xml, "var"); var_xml; var_xml = var_xml->next){
		varname = ezxml_attr(var_xml, "name");
		varpackages = ezxml_attr(var_xml, "packages");
		vararrgroup = ezxml_attr(var_xml, "array_group");
		varname_in_code = ezxml_attr(var_xml, "name_in_code");
		varname_in_output = ezxml_attr(var_xml, "name_in_output");
		skip_var = 0;

		if(!varname_in_code){
			varname_in_code = ezxml_attr(var_xml, "name");
		}

		for (var_xml2 = ezxml_child(var_arr_xml, "var"); var_xml2 && var_xml2 != var_xml; var_xml2 = var_xml2->next){
			// Check if the current array group has already been touched.
			vararrgroup2 = ezxml_attr(var_xml2, "array_group");

			if (strncmp(vararrgroup, vararrgroup2, 1024) == 0){
				skip_var = 1;
			}
		}

		if(!skip_var){

			fortprintf(fd, "! Starting group %s\n", vararrgroup);
			fortprintf(fd, "! Define constituent var %s\n", varname);
			fortprintf(fd, "! My Packages are %s\n", varpackages);

			// If no packages are defined, default to var_arr packages.
			if(varpackages == NULL){
				if(vararrpackages != NULL){
					varpackages = ezxml_attr(var_arr_xml, "packages");
				}
			}

			// Parse packages if they are defined
			sub_spacing[0] = '\0';
			if(varpackages){
				fortprintf(fd, "      if (");
				string = strdup(varpackages);
				tofree = string;
				token = strsep(&string, ";");
				fortprintf(fd, "%sActive", token);

				while( (token = strsep(&string, ";")) != NULL){
					fortprintf(fd, " .or. %sActive", token);
				}

				fortprintf(fd, ") then\n");
				snprintf(sub_spacing, 1024, "   ");
			}

			fortprintf(fd, "      %sindex_counter = index_counter + 1\n", sub_spacing);
			fortprintf(fd, "      %sif (associated(newSubPool)) then\n", sub_spacing);
			fortprintf(fd, "      %s   call mpas_pool_add_dimension(newSubPool, 'index_%s', index_counter)\n", sub_spacing, varname_in_code);
			fortprintf(fd, "      %send if\n", sub_spacing);
			fortprintf(fd, "      %sgroup_counter = group_counter + 1\n", sub_spacing);
			fortprintf(fd, "      %sif (.not. group_started) then\n", sub_spacing);
			fortprintf(fd, "      %s   group_start = index_counter\n", sub_spacing);
			fortprintf(fd, "      %s   if (associated(newSubPool)) then\n", sub_spacing);
			fortprintf(fd, "      %s      call mpas_pool_add_dimension(newSubPool, '%s_start', group_start)\n", sub_spacing, vararrgroup);
			fortprintf(fd, "      %s   end if\n", sub_spacing);
			fortprintf(fd, "      %s   group_started = .true.\n", sub_spacing);
			fortprintf(fd, "      %send if\n", sub_spacing);

			// If Packages are defined, write else clause
			if(varpackages){
				fortprintf(fd, "   %selse\n", sub_spacing);
				fortprintf(fd, "      %s  if (associated(newSubPool)) then\n", sub_spacing);
				fortprintf(fd, "      %s     call mpas_pool_add_dimension(newSubPool, 'index_%s', -1)\n", sub_spacing, varname_in_code);
				fortprintf(fd, "      %s  end if\n", sub_spacing);
				fortprintf(fd, "   %send if\n", sub_spacing);
			}

			// Add the rest of the variables from the current group.
			if(var_xml->next){
				for(var_xml2 = var_xml->next; var_xml2; var_xml2 = var_xml2->next){
					vararrgroup2 = ezxml_attr(var_xml2, "array_group");

					// var_xml2 is in the current array group
					if(strncmp(vararrgroup, vararrgroup2, 1024) == 0){
						varname = ezxml_attr(var_xml2, "name");
						varpackages = ezxml_attr(var_xml2, "packages");
						varname_in_code = ezxml_attr(var_xml2, "name_in_code");

						if(!varname_in_code){
							varname_in_code = ezxml_attr(var_xml2, "name");
						}


						// If no packages are defined, default to var_arr packages.
						if(varpackages == NULL){
							if(vararrpackages != NULL){
								varpackages = ezxml_attr(var_arr_xml, "packages");
							}
						}

						fortprintf(fd, "! Define constituent var %s\n", varname);
						fortprintf(fd, "! My packages are %s\n", varpackages);

						// Parse packages if they are defined
						sub_spacing[0] = '\0';
						if(varpackages){
							fortprintf(fd, "%sif (", spacing);
							string = strdup(varpackages);
							tofree = string;
							token = strsep(&string, ";");
							fortprintf(fd, "%sActive", token);

							while( (token = strsep(&string, ";")) != NULL){
								fortprintf(fd, " .or. %sActive", token);
							}

							fortprintf(fd, ") then\n");
							snprintf(sub_spacing, 1024, "   ");
						}

						fortprintf(fd, "      %sindex_counter = index_counter + 1\n", sub_spacing);
						fortprintf(fd, "      %sif (associated(newSubPool)) then\n", sub_spacing);
						fortprintf(fd, "      %s   call mpas_pool_add_dimension(newSubPool, 'index_%s', index_counter)\n", sub_spacing, varname_in_code);
						fortprintf(fd, "      %send if\n", sub_spacing);
						fortprintf(fd, "      %sgroup_counter = group_counter + 1\n", sub_spacing);
						fortprintf(fd, "      %sif (.not. group_started) then\n", sub_spacing);
						fortprintf(fd, "      %s   group_start = index_counter\n", sub_spacing);
						fortprintf(fd, "      %s   if (associated(newSubPool)) then\n", sub_spacing);
						fortprintf(fd, "      %s      call mpas_pool_add_dimension(newSubPool, '%s_start', group_start)\n", sub_spacing, vararrgroup);
						fortprintf(fd, "      %s   end if\n", sub_spacing);
						fortprintf(fd, "      %s   group_started = .true.\n", sub_spacing);
						fortprintf(fd, "      %send if\n", sub_spacing);

						// If Packages are defined, write else clause
						if(varpackages != NULL){
							fortprintf(fd, "   %selse\n", sub_spacing);
							fortprintf(fd, "   %s   if (associated(newSubPool)) then\n", sub_spacing);
							fortprintf(fd, "   %s      call mpas_pool_add_dimension(newSubPool, 'index_%s', -1)\n", sub_spacing, varname_in_code);
							fortprintf(fd, "   %s   end if\n", sub_spacing);
							fortprintf(fd, "   %send if\n", sub_spacing);
						}
					}
				}
			}

			fortprintf(fd, "      %sif (.not. group_started) then\n", sub_spacing);
			fortprintf(fd, "      %s   if (associated(newSubPool)) then\n", sub_spacing);
			fortprintf(fd, "      %s      call mpas_pool_add_dimension(newSubPool, '%s_start', -1)\n", sub_spacing, vararrgroup);
			fortprintf(fd, "      %s      call mpas_pool_add_dimension(newSubPool, '%s_end', -1)\n", sub_spacing, vararrgroup);
			fortprintf(fd, "      %s   end if\n", sub_spacing);
			fortprintf(fd, "      %selse\n", sub_spacing);
			fortprintf(fd, "      %s   group_started = .false.\n", sub_spacing);
			fortprintf(fd, "      %s   if (associated(newSubPool)) then\n", sub_spacing);
			fortprintf(fd, "      %s      call mpas_pool_add_dimension(newSubPool, '%s_end', index_counter)\n", sub_spacing, vararrgroup);
			fortprintf(fd, "      %s   end if\n", sub_spacing);
			fortprintf(fd, "      %send if\n", sub_spacing);
			fortprintf(fd, "! End of group       \n", vararrgroup);
		}
	}

	fortprintf(fd, "\n");

	// Setup constituent names
	fortprintf(fd, "      numConstituents = index_counter\n");
	fortprintf(fd, "      if (associated(newSubPool)) then\n");
	fortprintf(fd, "         call mpas_pool_add_dimension(newSubPool, 'num_%s', numConstituents)\n", vararrname);
	fortprintf(fd, "      end if\n");

	for(time_lev = 1; time_lev <= time_levs; time_lev++){
		if (time_levs > 1) {
			snprintf(pointer_name_arr, 1024, "%s(%d)", pointer_name, time_lev);
		} else {
			snprintf(pointer_name_arr, 1024, "%s", pointer_name);
		}
		fortprintf(fd, "! Defining time level %d\n", time_lev);
		fortprintf(fd, "      allocate( %s %% constituentNames(numConstituents) )\n", pointer_name_arr);
		fortprintf(fd, "      allocate( %s %% outputConstituentNames(numConstituents) )\n", pointer_name_arr);
		fortprintf(fd, "      %s %% fieldName = '%s'\n", pointer_name_arr , vararrname);
		fortprintf(fd, "      %s %% outputFieldName = '%s'\n", pointer_name_arr , vararrname);
		if (decomp != -1) {
			fortprintf(fd, "      %s %% isDecomposed = .true.\n", pointer_name_arr);
		} else {
			fortprintf(fd, "      %s %% isDecomposed = .false.\n", pointer_name_arr);
		}
		if (hasTime) {
			fortprintf(fd, "      %s %% hasTimeDimension = .true.\n", pointer_name_arr);
		} else {
			fortprintf(fd, "      %s %% hasTimeDimension = .false.\n", pointer_name_arr);
		}
		fortprintf(fd, "      %s %% isVarArray = .true.\n", pointer_name_arr);
		if(ndims > 0){
			if(persistence == SCRATCH){
				fortprintf(fd, "      %s %% isPersistent = .false.\n", pointer_name_arr);
				fortprintf(fd, "      %s %% isActive = .false.\n", pointer_name_arr);
			} else {
				fortprintf(fd, "      %s %% isPersistent = .true.\n", pointer_name_arr);
				fortprintf(fd, "      %s %% isActive = .false.\n", pointer_name_arr);
			}
		}
		fortprintf(fd, "\n");
		for(var_xml = ezxml_child(var_arr_xml, "var"); var_xml; var_xml = var_xml->next){
			varname = ezxml_attr(var_xml, "name");
			varname_in_code = ezxml_attr(var_xml, "name_in_code");

			if(!varname_in_code){
				varname_in_code = ezxml_attr(var_xml, "name");
			}

			varname_in_output = ezxml_attr(var_xml, "name_in_output");
			if(!varname_in_output){
				varname_in_output = ezxml_attr(var_xml, "name");
			}

			fortprintf(fd, "      if (associated(newSubPool)) then\n");
			fortprintf(fd, "         call mpas_pool_get_dimension(newSubPool, 'index_%s', const_index)\n", varname_in_code);
			fortprintf(fd, "      end if\n");
			fortprintf(fd, "      if (const_index > 0) then\n", spacing);
			fortprintf(fd, "         %s %% constituentNames(const_index) = '%s'\n", pointer_name_arr, varname);
			fortprintf(fd, "         %s %% outputConstituentNames(const_index) = '%s'\n", pointer_name_arr, varname_in_output);
			fortprintf(fd, "      end if\n", spacing);
		}

		fortprintf(fd, "\n");

		// Setup dimensions
		fortprintf(fd, "! Setup dimensions for       \n", vararrname);
		i = 1;
		fortprintf(fd, "      %s %% dimNames(%d) = 'num_%s'\n", pointer_name_arr, i, vararrname);

		string = strdup(vararrdims);
		tofree = string;
		token = strsep(&string, " ");

		if(strncmp(token, "Time", 1024) != 0){
			i++;
			if(strncmp(token, "nCells", 1024) == 0 || strncmp(token, "nEdges", 1024) == 0 || strncmp(token, "nVertices", 1024) == 0){
				fortprintf(fd, "      %s %% dimNames(%d) = '%s'\n", pointer_name_arr, i, token);
			} else {
				fortprintf(fd, "      %s %% dimNames(%d) = '%s'\n", pointer_name_arr, i, token);
			}
		}
		while( (token = strsep(&string, " ")) != NULL){
			if(strncmp(token, "Time", 1024) != 0){
				i++;
				if(strncmp(token, "nCells", 1024) == 0 || strncmp(token, "nEdges", 1024) == 0 || strncmp(token, "nVertices", 1024) == 0){
					fortprintf(fd, "      %s %% dimNames(%d) = '%s'\n", pointer_name_arr, i, token);
				} else {
					fortprintf(fd, "      %s %% dimNames(%d) = '%s'\n", pointer_name_arr, i, token);
				}
			}
		}
		free(tofree);

		fortprintf(fd, "\n");

		if ( ndims == 0 ) {
			fortprintf(fd, "      %s %% scalar = %s\n", pointer_name_arr, default_value);
		}
		fortprintf(fd, "      %s %% defaultValue = %s\n", pointer_name_arr, default_value);
		fortprintf(fd, "      allocate(%s %% attLists(size(%s %% constituentNames, dim=1)))\n", pointer_name_arr, pointer_name_arr);

		fortprintf(fd, "      do index_counter = 1, size(%s %% constituentNames, dim=1)\n", pointer_name_arr);
		fortprintf(fd, "         allocate(%s %% attLists(index_counter) %% attList)\n", pointer_name_arr);
		fortprintf(fd, "      end do\n");

		for(var_xml = ezxml_child(var_arr_xml, "var"); var_xml; var_xml = var_xml->next){
			varname = ezxml_attr(var_xml, "name");
			varname_in_code = ezxml_attr(var_xml, "name_in_code");
			vardesc = ezxml_attr(var_xml, "description");
			varunits = ezxml_attr(var_xml, "units");
			varbounds = ezxml_attr(var_xml, "bounds");
			varcellmeasures = ezxml_attr(var_xml, "cell_measures");
			varcellmethod = ezxml_attr(var_xml, "cell_method");
			varcoords = ezxml_attr(var_xml, "coordinates");
			varstdname = ezxml_attr(var_xml, "standard_name");

			if(!varname_in_code){
				varname_in_code = ezxml_attr(var_xml, "name");
			}

			fortprintf(fd, "      if (associated(newSubPool)) then\n");
			fortprintf(fd, "         call mpas_pool_get_dimension(newSubPool, 'index_%s', const_index)\n", varname_in_code);
			fortprintf(fd, "      end if\n");
			fortprintf(fd, "      if (const_index > 0) then\n", spacing);
			if ( vardesc != NULL ) {
				string = strdup(vardesc);
				tofree = string;

				token = strsep(&string, "'");
				sprintf(temp_str, "%s", token);

				while ( ( token = strsep(&string, "'") ) != NULL ) {
					sprintf(temp_str, "%s''%s", temp_str, token);
				}

				free(tofree);

				fortprintf(fd, "         call mpas_add_att(%s %% attLists(const_index) %% attList, 'long_name', '%s')\n", pointer_name_arr, temp_str);
			}

			if ( varunits != NULL ) {
				string = strdup(varunits);
				tofree = string;

				token = strsep(&string, "'");
				sprintf(temp_str, "%s", token);

				while ( ( token = strsep(&string, "'") ) != NULL ) {
					sprintf(temp_str, "%s''%s", temp_str, token);
				}

				free(tofree);

				fortprintf(fd, "         call mpas_add_att(%s %% attLists(const_index) %% attList, 'units', '%s')\n", pointer_name_arr, temp_str);
			}
		        if ( vararrmissing_value_mask != NULL && vararrmissingval != NULL ) {
				string = strdup(vararrmissing_value_mask);
				tofree = string;

				token = strsep(&string, "'");
				sprintf(temp_str, "%s", token);

				while ( ( token = strsep(&string, "'") ) != NULL ) {
					sprintf(temp_str, "%s''%s", temp_str, token);
				}

				free(tofree);

		                fortprintf(fd, "         %s %% maskName = '%s'\n", pointer_name_arr , temp_str);
				fortprintf(fd, "         call mpas_add_att(%s %% attLists(const_index) %% attList, 'missing_value_mask', '%s')\n", pointer_name_arr, temp_str);
                        } else {
                                fortprintf(fd, "         %s %% maskName = 'none'\n", pointer_name_arr , temp_str);
		        }

			if ( varbounds != NULL ) {
				string = strdup(varbounds);
				tofree = string;

				token = strsep(&string, "'");
				sprintf(temp_str, "%s", token);

				while ( ( token = strsep(&string, "'") ) != NULL ) {
					sprintf(temp_str, "%s''%s", temp_str, token);
				}

				free(tofree);

				fortprintf(fd, "         call mpas_add_att(%s %% attLists(const_index) %% attList, 'bounds', '%s')\n", pointer_name_arr, temp_str);
			} else if ( vararrbounds != NULL ) {
				string = strdup(vararrbounds);
				tofree = string;

				token = strsep(&string, "'");
				sprintf(temp_str, "%s", token);

				while ( ( token = strsep(&string, "'") ) != NULL ) {
					sprintf(temp_str, "%s''%s", temp_str, token);
				}

				free(tofree);

				fortprintf(fd, "         call mpas_add_att(%s %% attLists(const_index) %% attList, 'bounds', '%s')\n", pointer_name_arr, temp_str);
			}

			if ( varcellmeasures != NULL ) {
				string = strdup(varcellmeasures);
				tofree = string;

				token = strsep(&string, "'");
				sprintf(temp_str, "%s", token);

				while ( ( token = strsep(&string, "'") ) != NULL ) {
					sprintf(temp_str, "%s''%s", temp_str, token);
				}

				free(tofree);

				fortprintf(fd, "         call mpas_add_att(%s %% attLists(const_index) %% attList, 'cell_measures', '%s')\n", pointer_name_arr, temp_str);
			} else if ( vararrcellmeasures != NULL ) {
				string = strdup(vararrcellmeasures);
				tofree = string;

				token = strsep(&string, "'");
				sprintf(temp_str, "%s", token);

				while ( ( token = strsep(&string, "'") ) != NULL ) {
					sprintf(temp_str, "%s''%s", temp_str, token);
				}

				free(tofree);

				fortprintf(fd, "         call mpas_add_att(%s %% attLists(const_index) %% attList, 'cell_measures', '%s')\n", pointer_name_arr, temp_str);
			}

			if ( varcellmethod != NULL ) {
				string = strdup(varcellmethod);
				tofree = string;

				token = strsep(&string, "'");
				sprintf(temp_str, "%s", token);

				while ( ( token = strsep(&string, "'") ) != NULL ) {
					sprintf(temp_str, "%s''%s", temp_str, token);
				}

				free(tofree);

				fortprintf(fd, "         call mpas_add_att(%s %% attLists(const_index) %% attList, 'cell_method', '%s')\n", pointer_name_arr, temp_str);
			} else if ( vararrcellmethod != NULL ) {
				string = strdup(vararrcellmethod);
				tofree = string;

				token = strsep(&string, "'");
				sprintf(temp_str, "%s", token);

				while ( ( token = strsep(&string, "'") ) != NULL ) {
					sprintf(temp_str, "%s''%s", temp_str, token);
				}

				free(tofree);

				fortprintf(fd, "         call mpas_add_att(%s %% attLists(const_index) %% attList, 'cell_method', '%s')\n", pointer_name_arr, temp_str);
			}

			if ( varcoords != NULL ) {
				string = strdup(varcoords);
				tofree = string;

				token = strsep(&string, "'");
				sprintf(temp_str, "%s", token);

				while ( ( token = strsep(&string, "'") ) != NULL ) {
					sprintf(temp_str, "%s''%s", temp_str, token);
				}

				free(tofree);

				fortprintf(fd, "         call mpas_add_att(%s %% attLists(const_index) %% attList, 'coordinates', '%s')\n", pointer_name_arr, temp_str);
			} else if ( vararrcoords != NULL ) {
				string = strdup(vararrcoords);
				tofree = string;

				token = strsep(&string, "'");
				sprintf(temp_str, "%s", token);

				while ( ( token = strsep(&string, "'") ) != NULL ) {
					sprintf(temp_str, "%s''%s", temp_str, token);
				}

				free(tofree);

				fortprintf(fd, "         call mpas_add_att(%s %% attLists(const_index) %% attList, 'coordinates', '%s')\n", pointer_name_arr, temp_str);
			}

			if ( varstdname != NULL ) {
				string = strdup(varstdname);
				tofree = string;

				token = strsep(&string, "'");
				sprintf(temp_str, "%s", token);

				while ( ( token = strsep(&string, "'") ) != NULL ) {
					sprintf(temp_str, "%s''%s", temp_str, token);
				}

				free(tofree);

				fortprintf(fd, "         call mpas_add_att(%s %% attLists(const_index) %% attList, 'standard_name', '%s')\n", pointer_name_arr, temp_str);
			} else if ( vararrstdname != NULL ) {
				string = strdup(vararrstdname);
				tofree = string;

				token = strsep(&string, "'");
				sprintf(temp_str, "%s", token);

				while ( ( token = strsep(&string, "'") ) != NULL ) {
					sprintf(temp_str, "%s''%s", temp_str, token);
				}

				free(tofree);

				fortprintf(fd, "         call mpas_add_att(%s %% attLists(const_index) %% attList, 'standard_name', '%s')\n", pointer_name_arr, temp_str);
			}

			if ( vararrmissingval ) {
				fortprintf(fd, "         call mpas_add_att(%s %% attLists(const_index) %% attList, '_FillValue', %s)\n", pointer_name_arr, missing_value);
			}
			fortprintf(fd, "         %s %% missingValue = %s\n", pointer_name_arr, missing_value);
			fortprintf(fd, "         %s %% constituentNames(const_index) = '%s'\n", pointer_name_arr, varname);
			fortprintf(fd, "      end if\n", spacing);
		}


		fortprintf(fd, "      %s %% block => block\n", pointer_name_arr);
	}

	// Parse packages if they are defined
	fortprintf(fd, "\n");
	snprintf(spacing, 1024, "      ");
	if(!no_packages) {
		fortprintf(fd, "      if (");
		string = strdup(package_list);
		tofree = string;
		token = strsep(&string, ";");
		fortprintf(fd, "%sActive", token);
		sprintf(spacing, "   ");

		while( (token = strsep(&string, ";")) != NULL){
			fortprintf(fd, " .or. %sActive", token);
		}
		fortprintf(fd, ") then\n");
	}

	for(time_lev = 1; time_lev <= time_levs; time_lev++){
		if (time_levs > 1) {
			snprintf(pointer_name_arr, 1024, "%s(%d)", pointer_name, time_lev);
		} else {
			snprintf(pointer_name_arr, 1024, "%s", pointer_name);
		}
		fortprintf(fd, "      %s%s %% isActive = .true.\n", spacing, pointer_name_arr);
	}

	if (!no_packages) {
		fortprintf(fd, "      end if\n");
	}
	fortprintf(fd, "      call mpas_pool_add_field(newSubPool, '%s', %s)\n", vararrname_in_code, pointer_name);
	fortprintf(fd, "      call mpas_pool_add_field(block %% allFields, '%s', %s)\n", vararrname, pointer_name);
	fortprintf(fd, "\n");

	return 0;
}/*}}}*/


int parse_var(FILE *fd, ezxml_t registry, ezxml_t superStruct, ezxml_t currentVar, const char * corename)/*{{{*/
{
	ezxml_t struct_xml, var_xml, var_xml2;
	ezxml_t packages_xml, package_xml;
	ezxml_t streams_xml, stream_xml, streams_xml2, stream_xml2;

	const char *structtimelevs, *vartimelevs;
	const char *structname, *structlevs, *structpackages;
	const char *substructname;
	const char *varname, *varpersistence, *vartype, *vardims, *varunits, *vardesc, *vararrgroup, *varstreams, *vardefaultval, *varpackages, *varmissingval, *varmissing_value_mask;
	const char *varbounds, *varcellmeasures, *varcellmethod, *varcoords, *varstdname;
	const char *varname2, *vararrgroup2;
	const char *varname_in_code;
	const char *varname_in_output;
	const char *streamname, *streamname2;
	const char *packagename;

	int err;

	int iostreams;
	int i, skip_var, skip_stream;
	int time_lev, time_levs;
	int ndims, type, hasTime, decomp, in_stream;
	int persistence;
	char *string, *tofree, *token;
	char temp_str[1024];
	char pointer_name[1024];
	char pointer_name_arr[1024];
	char package_spacing[1024];
	char default_value[1024];
	char missing_value[1024];

	var_xml = currentVar;

	structname = ezxml_attr(superStruct, "name");
	structtimelevs = ezxml_attr(superStruct, "time_levs");

	// Define independent variables
	varname = ezxml_attr(var_xml, "name");
	vartype = ezxml_attr(var_xml, "type");
	vardims = ezxml_attr(var_xml, "dimensions");
	varpersistence = ezxml_attr(var_xml, "persistence");
	varpackages = ezxml_attr(var_xml, "packages");
	vardefaultval = ezxml_attr(var_xml, "default_value");
	vartimelevs = ezxml_attr(var_xml, "time_levs");
	varname_in_code = ezxml_attr(var_xml, "name_in_code");
	varunits = ezxml_attr(var_xml, "units");
	varbounds = ezxml_attr(var_xml, "bounds");
	varcellmeasures = ezxml_attr(var_xml, "cell_measures");
	varcellmethod = ezxml_attr(var_xml, "cell_method");
	varcoords = ezxml_attr(var_xml, "coordinates");
	varstdname = ezxml_attr(var_xml, "standard_name");
	vardesc = ezxml_attr(var_xml, "description");
	varmissingval = ezxml_attr(var_xml, "missing_value");
	varmissing_value_mask = ezxml_attr(var_xml, "missing_value_mask");

	if(!varname_in_code){
		varname_in_code = ezxml_attr(var_xml, "name");
	}

    varname_in_output = ezxml_attr(var_xml, "name_in_output");
    if(!varname_in_output){
        varname_in_output = ezxml_attr(var_xml, "name");
    }

	if(!vartimelevs){
		vartimelevs = ezxml_attr(superStruct, "time_levs");
	}

	if(vartimelevs){
		time_levs = atoi(vartimelevs);
		if(time_levs < 1){
			time_levs = 1;
		}
	} else {
		time_levs = 1;
	}

	persistence = check_persistence(varpersistence);

	fortprintf(fd, "! Define variable %s\n", varname);

	// Determine field type and default value.
	get_field_information(vartype, vardefaultval, default_value, varmissingval, missing_value, &type);

	// If a default_value is not specified, but a missing_value is, then set the
	// default_value (defaultValue) to missing_value
	if(!vardefaultval && varmissingval) {
		snprintf(default_value, 1024, "%s ! defaultValue taking specified missing_value", missing_value);
	}

	// Determine ndims, hasTime, and decomp type
	build_dimension_information(registry, var_xml, &ndims, &hasTime, &decomp);

	// Determine name of pointer for this field.
	set_pointer_name(type, ndims, pointer_name, time_levs);
	if (time_levs > 1) {
		fortprintf(fd, "      allocate(%s(%d))\n", pointer_name, time_levs);
	} else {
		fortprintf(fd, "      allocate(%s)\n", pointer_name);
	}

	for(time_lev = 1; time_lev <= time_levs; time_lev++){
		if (time_levs > 1) {
			snprintf(pointer_name_arr, 1024, "%s(%d)", pointer_name, time_lev);
		} else {
			snprintf(pointer_name_arr, 1024, "%s", pointer_name);
		}
		fortprintf(fd, "\n");
		fortprintf(fd, "! Setting up time level %d\n", time_lev);
		fortprintf(fd, "      %s %% fieldName = '%s'\n", pointer_name_arr, varname);
		fortprintf(fd, "      %s %% outputFieldName = '%s'\n", pointer_name_arr, varname_in_output);
		fortprintf(fd, "      %s %% isVarArray = .false.\n", pointer_name_arr);
		if (decomp != -1) {
			fortprintf(fd, "      %s %% isDecomposed = .true.\n", pointer_name_arr);
		} else {
			fortprintf(fd, "      %s %% isDecomposed = .false.\n", pointer_name_arr);
		}

		if(hasTime) {
			fortprintf(fd, "      %s %% hasTimeDimension = .true.\n", pointer_name_arr);
		} else {
			fortprintf(fd, "      %s %% hasTimeDimension = .false.\n", pointer_name_arr);
		}

		if(ndims > 0){
			if(persistence == SCRATCH){
				fortprintf(fd, "      %s %% isPersistent = .false.\n", pointer_name_arr);
				fortprintf(fd, "      %s %% isActive = .false.\n", pointer_name_arr);
			} else {
				fortprintf(fd, "      %s %% isPersistent = .true.\n", pointer_name_arr);
				fortprintf(fd, "      %s %% isActive = .false.\n", pointer_name_arr);
			}

			// Setup dimensions
			fortprintf(fd, "! Setting up dimensions\n");
			string = strdup(vardims);
			tofree = string;
			i = 1;
			token = strsep(&string, " ");
			if(strncmp(token, "Time", 1024) != 0){
				fortprintf(fd, "      %s %% dimNames(%d) = '%s'\n", pointer_name_arr, i, token);
				i++;
			}
			while( (token = strsep(&string, " ")) != NULL){
				if(strncmp(token, "Time", 1024) != 0){
					fortprintf(fd, "      %s %% dimNames(%d) = '%s'\n", pointer_name_arr, i, token);
					i++;
				}
			}
			free(tofree);
		}

		fortprintf(fd, "      %s %% defaultValue = %s\n", pointer_name_arr, default_value);
		if ( ndims == 0 ) {
			fortprintf(fd, "      %s %% scalar = %s\n", pointer_name_arr, default_value);
		}
		fortprintf(fd, "      allocate(%s %% attLists(1))\n", pointer_name_arr);
		fortprintf(fd, "      allocate(%s %% attLists(1) %% attList)\n", pointer_name_arr);

		if ( varunits != NULL ) {
			string = strdup(varunits);
			tofree = string;
			token = strsep(&string, "'");

			sprintf(temp_str, "%s", token);

			while ( ( token = strsep(&string, "'") ) != NULL ) {
				sprintf(temp_str, "%s''%s", temp_str, token);
			}

			free(tofree);

			fortprintf(fd, "      call mpas_add_att(%s %% attLists(1) %% attList, 'units', '%s')\n", pointer_name_arr, temp_str);
		}

		if ( varbounds != NULL ) {
			string = strdup(varbounds);
			tofree = string;
			token = strsep(&string, "'");

			sprintf(temp_str, "%s", token);

			while ( ( token = strsep(&string, "'") ) != NULL ) {
				sprintf(temp_str, "%s''%s", temp_str, token);
			}

			free(tofree);

			fortprintf(fd, "      call mpas_add_att(%s %% attLists(1) %% attList, 'bounds', '%s')\n", pointer_name_arr, temp_str);
		}

		if ( varcellmeasures != NULL ) {
			string = strdup(varcellmeasures);
			tofree = string;
			token = strsep(&string, "'");

			sprintf(temp_str, "%s", token);

			while ( ( token = strsep(&string, "'") ) != NULL ) {
				sprintf(temp_str, "%s''%s", temp_str, token);
			}

			free(tofree);

			fortprintf(fd, "      call mpas_add_att(%s %% attLists(1) %% attList, 'cell_measures', '%s')\n", pointer_name_arr, temp_str);
		}

		if ( varcellmethod != NULL ) {
			string = strdup(varcellmethod);
			tofree = string;
			token = strsep(&string, "'");

			sprintf(temp_str, "%s", token);

			while ( ( token = strsep(&string, "'") ) != NULL ) {
				sprintf(temp_str, "%s''%s", temp_str, token);
			}

			free(tofree);

			fortprintf(fd, "      call mpas_add_att(%s %% attLists(1) %% attList, 'cell_method', '%s')\n", pointer_name_arr, temp_str);
		}

		if ( varcoords != NULL ) {
			string = strdup(varcoords);
			tofree = string;
			token = strsep(&string, "'");

			sprintf(temp_str, "%s", token);

			while ( ( token = strsep(&string, "'") ) != NULL ) {
				sprintf(temp_str, "%s''%s", temp_str, token);
			}

			free(tofree);

			fortprintf(fd, "      call mpas_add_att(%s %% attLists(1) %% attList, 'coordinates', '%s')\n", pointer_name_arr, temp_str);
		}

		if ( varstdname != NULL ) {
			string = strdup(varstdname);
			tofree = string;
			token = strsep(&string, "'");

			sprintf(temp_str, "%s", token);

			while ( ( token = strsep(&string, "'") ) != NULL ) {
				sprintf(temp_str, "%s''%s", temp_str, token);
			}

			free(tofree);

			fortprintf(fd, "      call mpas_add_att(%s %% attLists(1) %% attList, 'standard_name', '%s')\n", pointer_name_arr, temp_str);
		}

		if ( vardesc != NULL ) {
			string = strdup(vardesc);
			tofree = string;
			token = strsep(&string, "'");

			sprintf(temp_str, "%s", token);

			while ( ( token = strsep(&string, "'") ) != NULL ) {
				sprintf(temp_str, "%s''%s", temp_str, token);
			}

			free(tofree);

			fortprintf(fd, "      call mpas_add_att(%s %% attLists(1) %% attList, 'long_name', '%s')\n", pointer_name_arr, temp_str);
		}

		if ( varmissing_value_mask != NULL && varmissingval != NULL ) {
			string = strdup(varmissing_value_mask);
			tofree = string;

			token = strsep(&string, "'");
			sprintf(temp_str, "%s", token);

			while ( ( token = strsep(&string, "'") ) != NULL ) {
				sprintf(temp_str, "%s''%s", temp_str, token);
			}

			free(tofree);

		        fortprintf(fd, "      %s %% maskName = '%s'\n", pointer_name_arr , temp_str);
			fortprintf(fd, "      call mpas_add_att(%s %% attLists(1) %% attList, 'missing_value_mask', '%s')\n", pointer_name_arr, temp_str);
                } else {
                        fortprintf(fd, "      %s %% maskName = 'none'\n", pointer_name_arr , temp_str);
		}

		if ( varmissingval != NULL ) {
			fortprintf(fd, "      call mpas_add_att(%s %% attLists(1) %% attList, '_FillValue', %s)\n", pointer_name_arr, missing_value);
		}
		fortprintf(fd, "      %s %% missingValue = %s\n", pointer_name_arr, missing_value);

		fortprintf(fd, "      %s %% block => block\n", pointer_name_arr);

	}

	// Parse packages if they are defined
	fortprintf(fd, "\n");
	package_spacing[0] = '\0';
	if(varpackages != NULL){
		fortprintf(fd, "      if (");
		string = strdup(varpackages);
		tofree = string;
		token = strsep(&string, ";");
		fortprintf(fd, "%sActive", token);
		sprintf(package_spacing, "   ");

		while( (token = strsep(&string, ";")) != NULL){
			fortprintf(fd, " .or. %sActive", token);
		}

		fortprintf(fd, ") then\n");
	}

	for(time_lev = 1; time_lev <= time_levs; time_lev++){
		if (time_levs > 1) {
			snprintf(pointer_name_arr, 1024, "%s(%d)", pointer_name, time_lev);
		} else {
			snprintf(pointer_name_arr, 1024, "%s", pointer_name);
		}
		fortprintf(fd, "      %s%s %% isActive = .true.\n", package_spacing, pointer_name_arr);
	}

	if(varpackages != NULL){
		fortprintf(fd, "      end if\n");
	}

	fortprintf(fd, "      call mpas_pool_add_field(newSubPool, '%s', %s)\n" , varname_in_code, pointer_name);
	fortprintf(fd, "      call mpas_pool_add_field(block %% allFields, '%s', %s)\n", varname, pointer_name);
	fortprintf(fd, "\n");

	return 0;
}/*}}}*/


int parse_struct(FILE *fd, ezxml_t registry, ezxml_t superStruct, int subpool, const char *parentname, const char * corename)/*{{{*/
{
	ezxml_t struct_xml, var_arr_xml, var_xml, var_xml2;
	ezxml_t struct_xml2;
	ezxml_t packages_xml, package_xml;

	const char *structname, *structlevs, *structpackages;
	const char *structname2;
	const char *substructname;
	const char *streamname, *streamname2;
	const char *packagename;
	const char *structnameincode;

	char *string, *tofree, *token;
	char spacing[1024];
	char core_string[1024];
	char pool_name[1024];
	char package_list[2048];

	int skip_struct, no_packages;
	int err;

	sprintf(core_string, "%s", corename);

	if(subpool){
		sprintf(pool_name, "%s_subpool", parentname);
	} else {
		sprintf(pool_name, "pool");
	}

	structname = ezxml_attr(superStruct, "name");
	structnameincode = ezxml_attr(superStruct, "name_in_code");

	if(!structnameincode){
		structnameincode = ezxml_attr(superStruct, "name");
	}

	structpackages = ezxml_attr(superStruct, "packages");

	// Extract all sub structs
	for (struct_xml = ezxml_child(superStruct, "var_struct"); struct_xml; struct_xml = struct_xml->next){
		err = parse_struct(fd, registry, struct_xml, 1, structname, corename);
	}

	fortprintf(fd, "   subroutine %s_generate_%s_%s(block, structPool, dimensionPool, packagePool)\n", core_string, pool_name, structname);
	fortprintf(fd, "      use mpas_derived_types\n");
	fortprintf(fd, "      use mpas_pool_routines\n");
	fortprintf(fd, "      use mpas_io_units\n");
	fortprintf(fd, "      use mpas_io, only : MPAS_REAL_FILLVAL, MPAS_INT_FILLVAL, MPAS_CHAR_FILLVAL\n");
	fortprintf(fd, "      implicit none\n");
	fortprintf(fd, "      type (block_type), intent(inout), pointer :: block\n");
	fortprintf(fd, "      type (mpas_pool_type), intent(inout) :: structPool\n");
	fortprintf(fd, "      type (mpas_pool_type), intent(inout) :: dimensionPool\n");
	fortprintf(fd, "      type (mpas_pool_type), intent(in) :: packagePool\n");
	write_field_pointer_arrays(fd);
	fortprintf(fd, "      type (mpas_pool_type), pointer :: newSubPool\n");
	fortprintf(fd, "      integer :: group_counter\n");
	fortprintf(fd, "      logical :: group_started\n");
	fortprintf(fd, "      integer :: group_start\n");
	fortprintf(fd, "      integer :: index_counter\n");
	fortprintf(fd, "      integer, pointer :: const_index\n");
	fortprintf(fd, "\n");

	// Need to define logicals for all packages
	for (packages_xml = ezxml_child(registry, "packages"); packages_xml; packages_xml = packages_xml->next){
		for (package_xml = ezxml_child(packages_xml, "package"); package_xml; package_xml = package_xml->next){
			packagename = ezxml_attr(package_xml, "name");

			fortprintf(fd, "      logical, pointer :: %sActive\n", packagename);
		}
	}

	fortprintf(fd, "\n");

	fortprintf(fd, "\n");
	fortprintf(fd, "      integer :: numConstituents\n");
	fortprintf(fd, "\n");

	fortprintf(fd, "      nullify(newSubPool)\n");

	fortprintf(fd, "      group_counter = -1\n");
	fortprintf(fd, "      group_started = .false.\n");
	fortprintf(fd, "      group_start = -1\n");

	// Need to get value of package flags
	for (packages_xml = ezxml_child(registry, "packages"); packages_xml; packages_xml = packages_xml->next){
		for (package_xml = ezxml_child(packages_xml, "package"); package_xml; package_xml = package_xml->next){
			packagename = ezxml_attr(package_xml, "name");

			fortprintf(fd, "      call mpas_pool_get_package(packagePool, '%sActive', %sActive)\n", packagename, packagename);
		}
	}

	fortprintf(fd, "\n");

	// Setup new pool to be added into structPool
	fortprintf(fd, "      allocate(newSubPool)\n");
	fortprintf(fd, "      call mpas_pool_create_pool(newSubPool)\n");
	fortprintf(fd, "      call mpas_pool_add_subpool(structPool, '%s', newSubPool)\n", structnameincode);
	fortprintf(fd, "      call mpas_pool_add_subpool(block %% allStructs, '%s', newSubPool)\n", structname);

	fortprintf(fd, "\n");

	// All sub-structs have been parsed and generated at this point. Time to generate this struct
	// Start by generating variable arrays
	for (var_arr_xml = ezxml_child(superStruct, "var_array"); var_arr_xml; var_arr_xml = var_arr_xml->next){
		parse_var_array(fd, registry, superStruct, var_arr_xml, corename);
	}


	// Define independent variables
	for (var_xml = ezxml_child(superStruct, "var"); var_xml; var_xml = var_xml->next){
		parse_var(fd, registry, superStruct, var_xml, corename);
	}

	fortprintf(fd, "\n");

	// Extract all sub structs
	for (struct_xml = ezxml_child(superStruct, "var_struct"); struct_xml; struct_xml = struct_xml->next){
		substructname = ezxml_attr(struct_xml, "name");
		fortprintf(fd, "      call %s_generate_%s_subpool_%s(block, newSubPool, dimensionPool, packagePool)\n", core_string, structname, substructname);
	}

	fortprintf(fd, "\n");
	fortprintf(fd, "      if (associated(newSubPool)) then\n");
	fortprintf(fd, "         call mpas_pool_add_config(newSubPool, 'on_a_sphere', block %% domain %% on_a_sphere)\n");
	fortprintf(fd, "         call mpas_pool_add_config(newSubPool, 'sphere_radius', block %% domain %% sphere_radius)\n");
	fortprintf(fd, "         call mpas_pool_add_config(newSubPool, 'is_periodic', block %% domain %% is_periodic)\n");
	fortprintf(fd, "         call mpas_pool_add_config(newSubPool, 'x_period', block %% domain %% x_period)\n");
	fortprintf(fd, "         call mpas_pool_add_config(newSubPool, 'y_period', block %% domain %% y_period)\n");
	fortprintf(fd, "      end if\n");
	fortprintf(fd, "\n");

	fortprintf(fd, "   end subroutine %s_generate_%s_%s\n", core_string, pool_name, structname);
	fortprintf(fd, "\n\n");

	return 0;
}/*}}}*/


int determine_struct_depth(int curLevel, ezxml_t superStruct){/*{{{*/
	ezxml_t subStruct;
	int max_depth, depth;

	max_depth = curLevel;

	for(subStruct = ezxml_child(superStruct, "var_struct"); subStruct; subStruct = subStruct->next){
		depth = determine_struct_depth(curLevel+1, subStruct);

		if(depth > max_depth){
			max_depth = depth;
		}
	}

	return max_depth;
}/*}}}*/


int generate_immutable_struct_contents(FILE *fd, const char *streamname, ezxml_t varstruct_xml){/*{{{*/

	ezxml_t var_xml, vararr_xml, substruct_xml;
	const char *optname, *optstream;

	/* Loop over fields looking for any that belong to the stream */
	for (vararr_xml = ezxml_child(varstruct_xml, "var"); vararr_xml; vararr_xml = vararr_xml->next) {
		optstream = ezxml_attr(vararr_xml, "streams");
		if (optstream != NULL && strstr(optstream, streamname) != NULL) {
			optname = ezxml_attr(vararr_xml, "name");
			fortprintf(fd, "   call MPAS_stream_mgr_add_field(manager, \'%s\', \'%s\', ierr=ierr)\n", streamname, optname);
		}
	}

	for (var_xml = ezxml_child(varstruct_xml, "var"); var_xml; var_xml = var_xml->next) {
		optstream = ezxml_attr(var_xml, "streams");
		if (optstream != NULL && strstr(optstream, streamname) != NULL) {
			optname = ezxml_attr(var_xml, "name");
			fortprintf(fd, "   call MPAS_stream_mgr_add_field(manager, \'%s\', \'%s\', ierr=ierr)\n", streamname, optname);
		}
	}

	for (substruct_xml = ezxml_child(varstruct_xml, "var_struct"); substruct_xml; substruct_xml = substruct_xml->next){
		generate_immutable_struct_contents(fd, streamname, substruct_xml);
	}

	return 0;
}/*}}}*/


/*********************************************************************************
 *
 *  Function: generate_immutable_streams
 *
 *  Generates the Fortran include file 'setup_immutable_streams.inc' that contains
 *  the subroutine mpas_generate_immutable_streams() responsible for making calls
 *  to the stream manager to define all immutable streams.
 *  The mpas_generate_immutable_streams() routine should be called after blocks
 *  have been allocated in the framework and after the stream manager has been
 *  initialized, but before any calls to generate mutable streams are made.
 *
 *********************************************************************************/
int generate_immutable_streams(ezxml_t registry){/*{{{*/
	ezxml_t streams_xml, stream_xml, var_xml, vararr_xml, varstruct_xml;
	ezxml_t substream_xml, matchstreams_xml, matchstream_xml;

	const char *optname, *opttype, *optvarname, *optstream, *optfilename, *optinterval_in, *optinterval_out, *optimmutable, *optpackages;
	const char *optstructname, *optsubstreamname, *optmatchstreamname, *optmatchimmutable;
	const char *corename;
	FILE *fd;

	char core_string[1024];

	corename = ezxml_attr(registry, "core_abbrev");

	sprintf(core_string, "%s", corename);

	fd = fopen("setup_immutable_streams.inc", "w+");

	fprintf(stderr, "---- GENERATING IMMUTABLE STREAMS ----\n");

	fortprintf(fd, "function %s_setup_immutable_streams(manager) result(iErr)\n\n", core_string);
	fortprintf(fd, "   use MPAS_derived_types, only : MPAS_streamManager_type, &\n");
	fortprintf(fd, "                                  MPAS_STREAM_INPUT_OUTPUT, MPAS_STREAM_INPUT, &\n");
	fortprintf(fd, "                                  MPAS_STREAM_OUTPUT, MPAS_STREAM_NONE, MPAS_STREAM_PROPERTY_IMMUTABLE\n");
	fortprintf(fd, "   use MPAS_stream_manager, only : MPAS_stream_mgr_create_stream, MPAS_stream_mgr_set_property, &\n");
	fortprintf(fd, "                                   MPAS_stream_mgr_add_field, MPAS_stream_mgr_add_pool\n");
	fortprintf(fd, "   use mpas_io_units\n\n");
	fortprintf(fd, "   implicit none\n\n");
	fortprintf(fd, "   type (MPAS_streamManager_type), pointer :: manager\n");
	fortprintf(fd, "   character (len=StrKIND) :: packages\n");
	fortprintf(fd, "   integer :: iErr\n\n");
	fortprintf(fd, "   iErr = 0\n\n");

	for (streams_xml = ezxml_child(registry, "streams"); streams_xml; streams_xml = streams_xml->next) {
		for (stream_xml = ezxml_child(streams_xml, "stream"); stream_xml; stream_xml = stream_xml->next) {

			optimmutable = ezxml_attr(stream_xml, "immutable");

			if (optimmutable != NULL && strcmp(optimmutable, "true") == 0) {

				optname = ezxml_attr(stream_xml, "name");
				opttype = ezxml_attr(stream_xml, "type");
				optfilename = ezxml_attr(stream_xml, "filename_template");

				/* create the stream */
				if (strstr(opttype, "input") != NULL && strstr(opttype, "output") != NULL)
					fortprintf(fd, "   call MPAS_stream_mgr_create_stream(manager, \'%s\', MPAS_STREAM_INPUT_OUTPUT, \'%s\', ierr=ierr)\n", optname, optfilename);
				else if (strstr(opttype, "input") != NULL)
					fortprintf(fd, "   call MPAS_stream_mgr_create_stream(manager, \'%s\', MPAS_STREAM_INPUT, \'%s\', ierr=ierr)\n", optname, optfilename);
				else if (strstr(opttype, "output") != NULL)
					fortprintf(fd, "   call MPAS_stream_mgr_create_stream(manager, \'%s\', MPAS_STREAM_OUTPUT, \'%s\', ierr=ierr)\n", optname, optfilename);
				else
					fortprintf(fd, "   call MPAS_stream_mgr_create_stream(manager, \'%s\', MPAS_STREAM_NONE, \'%s\', ierr=ierr)\n", optname, optfilename);


				/* Loop over streams listed within the stream (only use immutable streams) */
				for (substream_xml = ezxml_child(stream_xml, "stream"); substream_xml; substream_xml = ezxml_next(substream_xml)) {
					optsubstreamname = ezxml_attr(substream_xml, "name");
					optpackages = ezxml_attr(substream_xml, "packages");

					if (optpackages != NULL) {
						fortprintf(fd, "   write(packages,\'(a)\') \'%s\'\n", optpackages);
					}

					/* Find stream definition with matching name */
					for (matchstreams_xml = ezxml_child(registry, "streams"); matchstreams_xml; matchstreams_xml = matchstreams_xml->next){
						for (matchstream_xml = ezxml_child(matchstreams_xml, "stream"); matchstream_xml; matchstream_xml = matchstream_xml->next){
							optmatchstreamname = ezxml_attr(matchstream_xml, "name");
							optmatchimmutable = ezxml_attr(matchstream_xml, "immutable");

							if (optmatchstreamname != NULL && strcmp(optmatchstreamname, optsubstreamname) == 0){
								if (optmatchimmutable != NULL && strcmp(optmatchimmutable, "true") == 0) {
									/* Loop over fields listed within the stream */
									for (var_xml = ezxml_child(matchstream_xml, "var"); var_xml; var_xml = var_xml->next) {
										optvarname = ezxml_attr(var_xml, "name");
										if (optpackages != NULL)
											fortprintf(fd, "   call MPAS_stream_mgr_add_field(manager, \'%s\', \'%s\', packages=packages, ierr=ierr)\n", optname, optvarname);
										else
											fortprintf(fd, "   call MPAS_stream_mgr_add_field(manager, \'%s\', \'%s\', ierr=ierr)\n", optname, optvarname);

									}

									/* Loop over arrays of fields listed within the stream */
									for (vararr_xml = ezxml_child(matchstream_xml, "var_array"); vararr_xml; vararr_xml = vararr_xml->next) {
										optvarname = ezxml_attr(vararr_xml, "name");
										if (optpackages != NULL)
											fortprintf(fd, "   call MPAS_stream_mgr_add_field(manager, \'%s\', \'%s\', packages=packages, ierr=ierr)\n", optname, optvarname);
										else
											fortprintf(fd, "   call MPAS_stream_mgr_add_field(manager, \'%s\', \'%s\', ierr=ierr)\n", optname, optvarname);
									}

									/* Loop over var structs listed within the stream */
									for (varstruct_xml = ezxml_child(matchstream_xml, "var_struct"); varstruct_xml; varstruct_xml = varstruct_xml->next) {
										optstructname = ezxml_attr(varstruct_xml, "name");
										if (optpackages != NULL)
											fortprintf(fd, "   call MPAS_stream_mgr_add_pool(manager, \'%s\', \'%s\', packages=packages, ierr=ierr)\n", optname, optstructname);
										else
											fortprintf(fd, "   call MPAS_stream_mgr_add_pool(manager, \'%s\', \'%s\', ierr=ierr)\n", optname, optstructname);
									}

								} else {
									printf("ERROR: Immutable streams cannot contain mutable streams within them.\n");
									printf("ERROR:     Immutable stream \'%s\' contains a mutable stream \'%s\'.\n", optname, optsubstreamname);
									return 1;
								}
							}
						}
					}
				}

				/* Loop over var structs listed within the stream */
				for (varstruct_xml = ezxml_child(stream_xml, "var_struct"); varstruct_xml; varstruct_xml = ezxml_next(varstruct_xml)) {
					optstructname = ezxml_attr(varstruct_xml, "name");
					optpackages = ezxml_attr(varstruct_xml, "packages");

					if (optpackages != NULL) {
						fortprintf(fd, "   write(packages,\'(a)\') \'%s\'\n", optpackages);
						fortprintf(fd, "   call MPAS_stream_mgr_add_pool(manager, \'%s\', \'%s\', packages=packages, ierr=ierr)\n", optname, optstructname);
					}
					else {
						fortprintf(fd, "   call MPAS_stream_mgr_add_pool(manager, \'%s\', \'%s\', ierr=ierr)\n", optname, optstructname);
					}

				}


				/* Loop over arrays of fields listed within the stream */
				for (vararr_xml = ezxml_child(stream_xml, "var_array"); vararr_xml; vararr_xml = ezxml_next(vararr_xml)) {
					optvarname = ezxml_attr(vararr_xml, "name");
					optpackages = ezxml_attr(vararr_xml, "packages");

					if (optpackages != NULL) {
						fortprintf(fd, "   write(packages,\'(a)\') \'%s\'\n", optpackages);
						fortprintf(fd, "   call MPAS_stream_mgr_add_field(manager, \'%s\', \'%s\', packages=packages, ierr=ierr)\n", optname, optvarname);
					}
					else {
						fortprintf(fd, "   call MPAS_stream_mgr_add_field(manager, \'%s\', \'%s\', ierr=ierr)\n", optname, optvarname);
					}

				}

				/* Loop over fields listed within the stream */
				for (var_xml = ezxml_child(stream_xml, "var"); var_xml; var_xml = ezxml_next(var_xml)) {
					optvarname = ezxml_attr(var_xml, "name");
					optpackages = ezxml_attr(var_xml, "packages");

					if (optpackages != NULL) {
						fortprintf(fd, "   write(packages,\'(a)\') \'%s\'\n", optpackages);
						fortprintf(fd, "   call MPAS_stream_mgr_add_field(manager, \'%s\', \'%s\', packages=packages, ierr=ierr)\n", optname, optvarname);
					}
					else {
						fortprintf(fd, "   call MPAS_stream_mgr_add_field(manager, \'%s\', \'%s\', ierr=ierr)\n", optname, optvarname);
					}

				}

				/* Loop over fields looking for any that belong to the stream */
				for (varstruct_xml = ezxml_child(registry, "var_struct"); varstruct_xml; varstruct_xml = ezxml_next(varstruct_xml)) {
					generate_immutable_struct_contents(fd, optname, varstruct_xml);
				}

				fortprintf(fd, "   call MPAS_stream_mgr_set_property(manager, \'%s\', MPAS_STREAM_PROPERTY_IMMUTABLE, .true., ierr=ierr)\n\n", optname);
			}
		}
	}


	fortprintf(fd, "end function %s_setup_immutable_streams\n", core_string);

	fortprintf(fd, "\n\n");

	fclose(fd);

	return 0;
}/*}}}*/


int push_attributes(ezxml_t currentPosition){/*{{{*/
	ezxml_t child_xml, child_xml2;
	ezxml_t childStruct1, childStruct2, lastStruct;

	const char *name, *name2;
	const char *subname;
	const char *super_time_levs, *super_packages;
	const char *sub_time_levs, *sub_packages;

	int skip_struct;

	name = ezxml_attr(currentPosition, "name");

	// Iterate over var_arrays
	for(child_xml = ezxml_child(currentPosition, "var_array"); child_xml; child_xml = child_xml->next){
		super_time_levs = ezxml_attr(currentPosition, "time_levs");
		super_packages = ezxml_attr(currentPosition, "packages");
		subname = ezxml_attr(child_xml, "name");
		sub_time_levs = ezxml_attr(child_xml, "time_levs");
		sub_packages = ezxml_attr(child_xml, "packages");

		if(!sub_time_levs && super_time_levs){
			child_xml = ezxml_set_attr(child_xml, "time_levs", super_time_levs);
		}

		if(!sub_packages && super_packages){
			child_xml = ezxml_set_attr(child_xml, "packages", super_packages);
		}

		// Iterate over vars in var_array
		for(child_xml2 = ezxml_child(child_xml, "var"); child_xml2; child_xml2 = child_xml2->next){
			super_packages = ezxml_attr(child_xml, "packages");
			sub_packages = ezxml_attr(child_xml2, "packages");

			if(!sub_packages && super_packages){
				child_xml2 = ezxml_set_attr(child_xml2, "packages", super_packages);
			}
		}
	}

	// Iterate over vars
	for(child_xml = ezxml_child(currentPosition, "var"); child_xml; child_xml = child_xml->next){
		super_packages = ezxml_attr(currentPosition, "packages");
		super_time_levs = ezxml_attr(currentPosition, "time_levs");
		subname = ezxml_attr(child_xml, "name");
		sub_time_levs = ezxml_attr(child_xml, "time_levs");
		sub_packages = ezxml_attr(child_xml, "packages");

		if(!sub_time_levs && super_time_levs){
			child_xml = ezxml_set_attr(child_xml, "time_levs", super_time_levs);
		}

		if(!sub_packages && super_packages){
			child_xml = ezxml_set_attr(child_xml, "packages", super_packages);
		}
	}

	// Iterate over var structs
	for(child_xml = ezxml_child(currentPosition, "var_struct"); child_xml; child_xml = child_xml->next){
		super_packages = ezxml_attr(currentPosition, "packages");
		super_time_levs = ezxml_attr(currentPosition, "time_levs");

		subname = ezxml_attr(child_xml, "name");
		sub_time_levs = ezxml_attr(child_xml, "time_levs");
		sub_packages = ezxml_attr(child_xml, "packages");

		if(!sub_time_levs && super_time_levs){
			child_xml = ezxml_set_attr(child_xml, "time_levs", super_time_levs);
		}

		if(!sub_packages && super_packages){
			child_xml = ezxml_set_attr(child_xml, "packages", super_packages);
		}

		push_attributes(child_xml);
	}

	return 0;
}/*}}}*/


int merge_structs_and_var_arrays(ezxml_t currentPosition){/*{{{*/
	ezxml_t old_child, new_child;
	ezxml_t next_child;
	ezxml_t childStruct1, childStruct2, lastStruct;

	const char *name, *name2;
	const char *subname;

	int skip_struct;

	// Merge var_structs
	for(childStruct1 = ezxml_child(currentPosition, "var_struct"); childStruct1; childStruct1 = childStruct1->next){
		name = ezxml_attr(childStruct1, "name");

		skip_struct = 0;
		for(childStruct2 = ezxml_child(currentPosition, "var_struct"); childStruct2 != childStruct1 && childStruct2; childStruct2 = childStruct2->next){
			name2 = ezxml_attr(childStruct2, "name");

			if(strcmp(name, name2) == 0){
				skip_struct = 1;
			}
		}

		if(!skip_struct && childStruct1->next){
			lastStruct = childStruct1;
			for(childStruct2 = childStruct1->next; childStruct2; childStruct2 = childStruct2->next){
				name2 = ezxml_attr(childStruct2, "name");

				if(strcmp(name, name2) == 0){
					// Merge children into childStruct1, and "remove" childStruct2
					for(old_child = ezxml_child(childStruct2, "var"); old_child; old_child = next_child){
						next_child = ezxml_next(old_child);
						new_child = ezxml_insert(old_child, childStruct1, strlen(childStruct1->txt));
					}

					for(old_child = ezxml_child(childStruct2, "var_array"); old_child; old_child = next_child){
						next_child = ezxml_next(old_child);
						new_child = ezxml_insert(old_child, childStruct1, strlen(childStruct1->txt));
					}

					for(old_child = ezxml_child(childStruct2, "var_struct"); old_child; old_child = next_child){
						next_child = ezxml_next(old_child);
						new_child = ezxml_insert(old_child, childStruct1, strlen(childStruct1->txt));
					}

					// Remove childStruct2
					lastStruct->next = childStruct2->next;
					free(childStruct2);
					childStruct2 = lastStruct;
				} else {
					lastStruct = childStruct2;
				}
			}
		}
	}

	// Merge var_arrays
	for(childStruct1 = ezxml_child(currentPosition, "var_array"); childStruct1; childStruct1 = childStruct1->next){
		name = ezxml_attr(childStruct1, "name");

		skip_struct = 0;
		for(childStruct2 = ezxml_child(currentPosition, "var_array"); childStruct2 && childStruct2 != childStruct1; childStruct2 = childStruct2->next){
			name2 = ezxml_attr(childStruct2, "name");

			if(strcmp(name, name2) == 0){
				skip_struct = 1;
			}
		}

		if(!skip_struct && childStruct1->next){
			lastStruct = childStruct1;
			for(childStruct2 = childStruct1->next; childStruct2; childStruct2 = childStruct2->next){
				name2 = ezxml_attr(childStruct2, "name");

				if(strcmp(name, name2) == 0){
					// Merge var_arrays and remove childStruct2
					for(old_child = ezxml_child(childStruct2, "var"); old_child; old_child = next_child){
						next_child = ezxml_next(old_child);
						new_child = ezxml_insert(old_child, childStruct1, strlen(childStruct1->txt));
					}

					lastStruct->next = childStruct2->next;
					free(childStruct2);
					childStruct2 = lastStruct;
				} else {
					lastStruct = childStruct2;
				}
			}
		}
	}

	for(childStruct1 = ezxml_child(currentPosition, "var_struct"); childStruct1; childStruct1 = childStruct1->next){
		merge_structs_and_var_arrays(childStruct1);
	}

	return 0;
}/*}}}*/


int merge_streams(ezxml_t registry){/*{{{*/
	ezxml_t old_child, new_child, tmp_child;
	ezxml_t childStream1, childStream2, lastStream;
	ezxml_t includeStream;

	ezxml_t streamsBlock, streamsBlock2;

	const char *name, *name2;
	const char *subname;

	int skip_stream;

	// First, merge all streams blocks. Regardless of nested stream names.
	streamsBlock = ezxml_child(registry, "streams");
	while(streamsBlock->next){
		for(childStream1 = ezxml_child(streamsBlock->next, "stream"); childStream1; childStream1){
			if (childStream1->next){
				lastStream = childStream1->next;
			} else {
				lastStream = NULL;
			}
			name = ezxml_attr(childStream1, "name");
			new_child = ezxml_insert(childStream1, streamsBlock, strlen(streamsBlock->txt));

			childStream1 = lastStream;

		}

		lastStream = streamsBlock->next;
		streamsBlock->next = streamsBlock->next->next;
		free(lastStream);
	}


	// Now, merge all streams with the same name, within streamsBlock
	streamsBlock = ezxml_child(registry, "streams");
	for(childStream1 = ezxml_child(streamsBlock, "stream"); childStream1; childStream1 = childStream1->next){
		name = ezxml_attr(childStream1, "name");

		skip_stream = 0;

		for(childStream2 = ezxml_child(streamsBlock, "stream"); childStream2 && childStream2 != childStream1; childStream2 = childStream2->next){
			name2 = ezxml_attr(childStream2, "name");

			if(strcmp(name, name2) == 0){
				skip_stream = 1;
			}
		}

		if(!skip_stream && childStream1->next){
			lastStream = childStream1;
			for(childStream2 = childStream1->next; childStream2; childStream2 = childStream2->next){
				name2 = ezxml_attr(childStream2, "name");

				if(strcmp(name, name2) == 0){
					// Merge child vars
					for(old_child = ezxml_child(childStream2, "var"); old_child; old_child){
						if(old_child->next){
							tmp_child = old_child->next;
						} else {
							tmp_child = NULL;
						}
						new_child = ezxml_insert(old_child, childStream1, strlen(childStream1->txt));

						old_child = tmp_child;
					}

					for(old_child = ezxml_child(childStream2, "var_array"); old_child; old_child){
						if(old_child->next){
							tmp_child = old_child->next;
						} else {
							tmp_child = NULL;
						}
						new_child = ezxml_insert(old_child, childStream1, strlen(childStream1->txt));

						old_child = tmp_child;
					}

					for(old_child = ezxml_child(childStream2, "var_struct"); old_child; old_child){
						if(old_child->next){
							tmp_child = old_child->next;
						} else {
							tmp_child = NULL;
						}
						new_child = ezxml_insert(old_child, childStream1, strlen(childStream1->txt));

						old_child = tmp_child;
					}

					for(old_child = ezxml_child(childStream2, "stream"); old_child; old_child){
						if(old_child->next){
							tmp_child = old_child->next;
						} else {
							tmp_child = NULL;
						}
						new_child = ezxml_insert(old_child, childStream1, strlen(childStream1->txt));

						old_child = tmp_child;
					}




					lastStream->next = childStream2->next;
					free(childStream2);
					childStream2 = lastStream;
				} else {
					lastStream = childStream2;
				}
			}
		}
	}

	return 0;
}/*}}}*/


int parse_structs_from_registry(ezxml_t registry)/*{{{*/
{
	ezxml_t structs_xml, var_arr_xml, var_xml;
	ezxml_t packages_xml, package_xml;

	const char *corename, *packagename, *structname, *structpackages;
	FILE *fd;
	int err;

	char core_string[1024];
	char spacing[1024];
	char package_list[2048];

	int no_packages;

	char *string, *tofree, *token;

	corename = ezxml_attr(registry, "core_abbrev");

	sprintf(core_string, "%s", corename);

	fd = fopen("structs_and_variables.inc", "w+");

	for (structs_xml = ezxml_child(registry, "var_struct"); structs_xml; structs_xml = structs_xml->next){
		err = parse_struct(fd, registry, structs_xml, 0, "\0", corename);
	}

	fortprintf(fd, "   subroutine %s_generate_structs(block, structPool, dimensionPool, packagePool)\n", core_string);
	fortprintf(fd, "      use mpas_derived_types\n");
	fortprintf(fd, "      use mpas_io_units\n");
	fortprintf(fd, "      implicit none\n");
	fortprintf(fd, "      type (block_type), pointer, intent(inout) :: block\n");
	fortprintf(fd, "      type (mpas_pool_type), intent(inout) :: structPool\n");
	fortprintf(fd, "      type (mpas_pool_type), intent(inout) :: dimensionPool\n");
	fortprintf(fd, "      type (mpas_pool_type), intent(in) :: packagePool\n");

	fortprintf(fd, "\n");

	for (structs_xml = ezxml_child(registry, "var_struct"); structs_xml; structs_xml = structs_xml->next){
		structname = ezxml_attr(structs_xml, "name");

		fortprintf(fd, "      call %s_generate_pool_%s(block, structPool, dimensionPool, packagePool)\n", core_string, structname);

		fortprintf(fd, "\n");
	}


	fortprintf(fd, "   end subroutine %s_generate_structs\n", core_string);

	fclose(fd);

	return 0;
}/*}}}*/


