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
#include "registry_types.h"
#include "gen_inc.h"
#include "ezxml/ezxml.h"

int getword(FILE *, char *);
int is_integer_constant(char *);
void sort_vars(struct variable *);
void sort_group_vars(struct group_list *);
int parse_reg_xml(ezxml_t registry, struct namelist **nls, struct dimension ** dims, struct variable ** vars, struct group_list ** groups, struct package ** pkgs, char * modelname, char * corename, char * version);
int validate_reg_xml(ezxml_t registry);
char * check_packages(ezxml_t registry, char * dims);
char * check_dimensions(ezxml_t registry, char * dims);
char * check_streams(char * streams);
int check_persistence(const char * persistence);

int main(int argc, char ** argv)/*{{{*/
{
	FILE * regfile;
	struct namelist * nls;
	struct dimension * dims;
	struct variable * vars;
	struct group_list * groups;
	struct package * pkgs;

	char *modelname, *corename, *version;

	modelname = (char *)malloc(sizeof(char)*1024);
	corename = (char *)malloc(sizeof(char)*1024);
	version = (char *)malloc(sizeof(char)*1024);

	if (argc != 2) {
		fprintf(stderr,"Reading registry file from standard input\n");
		regfile = stdin;
	}
	else if (!(regfile = fopen(argv[1], "r"))) {
		fprintf(stderr,"\nError: Could not open file %s for reading.\n\n", argv[1]);
		return 1;
	}   

	nls = NULL;
	dims = NULL;
	vars = NULL;

	ezxml_t registry = ezxml_parse_fp(regfile);

	if (validate_reg_xml(registry)) {
		fprintf(stderr, "Validation failed.....\n");
		return 1;
	}

	if (parse_reg_xml(registry, &nls, &dims, &vars, &groups, &pkgs, modelname, corename, version)) {
		fprintf(stderr, "Parsing failed.....\n");
		return 1;
	}

	sort_vars(vars);
	sort_group_vars(groups);

	write_default_namelist(nls, corename);
	gen_history_attributes(modelname, corename, version);
	gen_namelists(nls);
	gen_field_defs(groups, vars, dims);
	gen_reads(groups, vars, dims);
	gen_writes(groups, vars, dims, nls);
	gen_packages(pkgs);

	free(modelname);
	free(corename);
	free(version);

	return 0;
}/*}}}*/

int validate_reg_xml(ezxml_t registry)/*{{{*/
{
	ezxml_t dims_xml, dim_xml;
	ezxml_t structs_xml, var_arr_xml, var_xml;
	ezxml_t nmlrecs_xml, nmlopt_xml;
	ezxml_t streams_xml, stream_xml;

	const char *dimname, *dimunits, *dimdesc, *dimdef;
	const char *nmlrecname, *nmlrecindef;
	const char *nmloptname, *nmlopttype, *nmloptval, *nmloptunits, *nmloptdesc, *nmloptposvals, *nmloptindef;
	const char *structname, *structlevs, *structpackages;
	const char *vararrname, *vararrtype, *vararrdims, *vararrpersistence, *vararrpackages;
	const char *varname, *varpersistence, *vartype, *vardims, *varunits, *vardesc, *vararrgroup, *varstreams, *varpackages;
	const char *varname_in_code;
	const char *const_model, *const_core, *const_version;
	const char *streamname;
	const char *streamtype;

	char *string, *err_string;
	char name_holder[1024];

	int found;
	int persistence;

	// Get model information
	const_model = ezxml_attr(registry, "model");
	const_core = ezxml_attr(registry, "core");
	const_version = ezxml_attr(registry, "version");

	if(const_model == NULL)
		fprintf(stderr,"Warning: Model attribute missing in registry declaration.\n");

	if(const_core == NULL)
		fprintf(stderr,"Warning: Core attribute missing in registry declaration.\n");

	if(const_version == NULL)
		fprintf(stderr,"Warning: Version attribute missing in registry declaration.\n");

	// Validate default streams
	for (streams_xml = ezxml_child(registry, "streams"); streams_xml; streams_xml = streams_xml->next) {
		for (stream_xml = ezxml_child(streams_xml, "stream"); stream_xml; stream_xml = stream_xml->next) {
			streamname = ezxml_attr(stream_xml, "name");
			streamtype = ezxml_attr(stream_xml, "type");
			if (streamname == NULL) {
				fprintf(stderr,"ERROR: Stream specification missing \"name\" attribute.\n");
				return 1;
			}
			else if (streamtype == NULL) {
				fprintf(stderr,"ERROR: Stream specification missing \"type\" attribute.\n");
				return 1;
			}
			else {
				for (var_xml = ezxml_child(stream_xml, "var"); var_xml; var_xml = var_xml->next) {
					varname = ezxml_attr(var_xml, "name");
					if (varname == NULL) {
						fprintf(stderr,"ERROR: Variable field in stream \"%s\" specification missing \"name\" attribute.\n", streamname);
						return 1;
					}
				}
			}
		}
	}

	// Validate Namelist Records
	for (nmlrecs_xml = ezxml_child(registry, "nml_record"); nmlrecs_xml; nmlrecs_xml = nmlrecs_xml->next){
		nmlrecname = ezxml_attr(nmlrecs_xml, "name");

		if (nmlrecname == NULL){
			fprintf(stderr,"ERROR: Namelist record missing name attribute.\n");
			return 1;
		}

		for (nmlopt_xml = ezxml_child(nmlrecs_xml, "nml_option"); nmlopt_xml; nmlopt_xml = nmlopt_xml->next){
			nmloptname = ezxml_attr(nmlopt_xml, "name");
			nmlopttype = ezxml_attr(nmlopt_xml, "type");
			nmloptval = ezxml_attr(nmlopt_xml, "default_value");
			nmloptunits = ezxml_attr(nmlopt_xml, "units");
			nmloptdesc = ezxml_attr(nmlopt_xml, "description");
			nmloptposvals = ezxml_attr(nmlopt_xml, "possible_values");

			if (nmloptname == NULL){
				fprintf(stderr,"ERROR: Namelist option missing name attribute in record %s.\n", nmlrecname);
				return 1;
			}

			if (nmlopttype == NULL){
				fprintf(stderr,"ERROR: Namelist option %s missing type attribute.\n", nmloptname);
				return 1;
			} else if (strcasecmp("logical", nmlopttype) != 0 && strcasecmp("real", nmlopttype) != 0 &&
					strcasecmp("integer", nmlopttype) != 0 && strcasecmp("character", nmlopttype) != 0) {
				fprintf(stderr,"ERROR: Type of namelist option %s doesn't equal one of logical, real, character, or integer.\n", nmloptname);
				return 1;
			}

			if (nmloptval == NULL){
				fprintf(stderr,"ERROR: Default value missing for namelist option %s.\n", nmloptname);
				return 1;
			}
		}
	}

	// Validate Dimensions
	for (dims_xml = ezxml_child(registry, "dims"); dims_xml; dims_xml = dims_xml->next){
		for (dim_xml = ezxml_child(dims_xml, "dim"); dim_xml; dim_xml = dim_xml->next){
			dimname = ezxml_attr(dim_xml, "name");
			dimdef = ezxml_attr(dim_xml, "definition");	
			dimunits = ezxml_attr(dim_xml, "units");
			dimdesc = ezxml_attr(dim_xml, "description");

			if (dimname == NULL){
				fprintf(stderr,"ERROR: Name missing for dimension.\n");
				return 1;
			}

			if (dimdef != NULL){
				if (strncmp(dimdef, "namelist:", 9) == 0){
					found = 0;
					snprintf(name_holder, 1024, "%s",dimdef);
					snprintf(name_holder, 1024, "%s",(name_holder)+9);
					for (nmlrecs_xml = ezxml_child(registry, "nml_record"); nmlrecs_xml; nmlrecs_xml = nmlrecs_xml->next){
						nmlrecindef = ezxml_attr(nmlrecs_xml, "in_defaults");

						if(nmlrecindef != NULL){
							if(strncmp(nmlrecindef, "true", 1024) != 0 && strncmp(nmlrecindef, "false", 1024) != 0){
								fprintf(stderr, "ERROR: Namelist record %s has an invalid value for in_defaults attribute. Valide values are true or false.\n", nmlrecname);
							}
						}
						for (nmlopt_xml = ezxml_child(nmlrecs_xml, "nml_option"); nmlopt_xml; nmlopt_xml = nmlopt_xml->next){
							nmloptname = ezxml_attr(nmlopt_xml, "name");
							nmlopttype = ezxml_attr(nmlopt_xml, "type");
							nmloptindef = ezxml_attr(nmlopt_xml, "in_defaults");


							if(nmloptindef != NULL){
								if(strncmp(nmloptindef, "true", 1024) != 0 && strncmp(nmloptindef, "false", 1024) != 0){
									fprintf(stderr, "ERROR: Namelist option %s in record %s has an invalid value for in_defaults attribute. Valide values are true or false.\n", nmloptname, nmlrecname);
								}
							}

							if (strncmp(name_holder, nmloptname, 1024) == 0){
								if (strcasecmp("integer", nmlopttype) != 0){
									fprintf(stderr, "ERROR: Namelist variable %s must be an integer for namelist-derived dimension %s.\n", nmloptname, dimname);
									return 1;
								}

								found = 1;
							}
						}
					}

					if (!found){
						fprintf(stderr, "ERROR: Namelist variable %s not found for namelist-derived dimension %s\n", name_holder, dimname);
						return 1;
					}
				}
			}
		}
	}

	// Validate Variable Structures
	for(structs_xml = ezxml_child(registry, "var_struct"); structs_xml; structs_xml = structs_xml->next){
		structname = ezxml_attr(structs_xml, "name");
		structlevs = ezxml_attr(structs_xml, "time_levs");
		structpackages = ezxml_attr(structs_xml, "packages");

		if (structname == NULL){
			fprintf(stderr,"ERROR: Name missing for var_struct.\n");
			return 1;
		}

		if (structlevs == NULL){
			fprintf(stderr,"ERROR: time_levs attribute missing for var_struct %s.\n", structname);
			return 1;
		}

		if (structpackages != NULL) {
			string = strdup(structpackages);
			err_string = check_packages(registry, string);
			free(string);

			if (err_string != NULL){
				fprintf(stderr, "ERROR: Package %s used on var_struct %s is not defined.\n", err_string, structname);
				return 1;
			}
		}

		// Validate variable arrays
		for(var_arr_xml = ezxml_child(structs_xml, "var_array"); var_arr_xml; var_arr_xml = var_arr_xml->next){
			vararrname = ezxml_attr(var_arr_xml, "name");
			vararrtype = ezxml_attr(var_arr_xml, "type");
			vararrdims = ezxml_attr(var_arr_xml, "dimensions");
			vararrpersistence = ezxml_attr(var_arr_xml, "persistence");
			vararrpackages = ezxml_attr(var_arr_xml, "packages");

			if (vararrname == NULL){
				fprintf(stderr,"ERROR: Name attribute missing for var_array in var_struct %s.\n", structname);
				return 1;
			}

			if (vararrtype == NULL){
				fprintf(stderr,"ERROR: Type attribute missing for var_array %s in var_struct %s.\n", vararrname, structname);
				return 1;
			} else if (strcasecmp("logical", vararrtype) != 0 && strcasecmp("real", vararrtype) != 0 &&
					strcasecmp("integer", vararrtype) != 0 && strcasecmp("text", vararrtype) != 0) {
				fprintf(stderr,"ERROR: Type attribute on var_array %s in var_struct %s is not equal to one of logical, real, integer, or text.\n", vararrname, structname);
				return 1;
			}

			if (vararrdims == NULL){
				fprintf(stderr,"ERROR: Dimensions attribute missing for var_array %s in var_struct %s.\n", vararrname, structname);
				return 1;
			} else { 
				string = strdup(vararrdims);
				err_string = check_dimensions(registry, string);
				free(string);

				if (err_string != NULL){
					fprintf(stderr,"ERROR: Dimension %s on var_array %s in var_struct %s is not defined.\n", err_string, vararrname, structname);
					return 1;
				}
			}

			persistence = PERSISTENT;
			if (vararrpersistence != NULL){
				persistence = check_persistence(vararrpersistence);

				if(persistence == -1) {
					fprintf(stderr, "\ton var_array %s in var_struct %s.\n", vararrname, structname);
					return -1;
				}
			}

			if(persistence == SCRATCH && vararrpackages != NULL){
				fprintf(stderr, "ERROR: Packages attribute not allowed on scratch var_array %s in var_struct %s.\n", vararrname, structname);
				return -1;
			} else if (persistence == SCRATCH && vararrpackages == NULL && structpackages != NULL) {
				fprintf(stderr, "ERROR: Packages attribute inherited from var_struct %s not allowed on scratch var_array %s in var_struct %s.\n", structname, vararrname, structname);
				return -1;
			} else if (persistence == PERSISTENT && vararrpackages != NULL){
				string = strdup(vararrpackages);
				err_string = check_packages(registry, string);
				free(string);

				if (err_string != NULL){
					fprintf(stderr, "ERROR: Package %s used on var_array %s in var_struct %s is not defined.\n", err_string, vararrname, structname);
					return 1;
				}
			}


			// Validate variables in variable arrays
			for(var_xml = ezxml_child(var_arr_xml, "var"); var_xml; var_xml = var_xml->next){
				varname = ezxml_attr(var_xml, "name");
				varunits = ezxml_attr(var_xml, "units");
				vardesc = ezxml_attr(var_xml, "description");
				vararrgroup = ezxml_attr(var_xml, "array_group");
				varname_in_code = ezxml_attr(var_xml, "name_in_code");
				varpackages = ezxml_attr(var_xml, "packages");

				if (varname == NULL) {
					fprintf(stderr,"ERROR: Name missing for constituent variable in var_array %s in var_struct %s.\n", vararrname, structname);
					return 1;
				}

				if (vararrgroup == NULL){
					fprintf(stderr,"ERROR: Array group attribute missing for constituent variable %s in var_array %s in var_struct %s.\n", varname, vararrname, structname);
					return 1;
				}

				if (persistence == SCRATCH && vararrpackages != NULL) {
					fprintf(stderr, "ERROR: Packages attribute not allowed on constituent variable %s within scratch var_srray %s in var_struct %s.\n", varname, vararrname, structname);
					return 1;
				}

				if(varpackages != NULL){
					string = strdup(varpackages);
					err_string = check_packages(registry, string);
					free(string);

					if (err_string != NULL){
						fprintf(stderr, "ERROR: Package %s used on constituent variable %s in var_array %s var_struct %s is not defined.\n", err_string, varname, vararrname, structname);
						return 1;
					}
				}

			}
		}

		for(var_xml = ezxml_child(structs_xml, "var"); var_xml; var_xml = var_xml->next){
			varname = ezxml_attr(var_xml, "name");
			varpersistence = ezxml_attr(var_xml, "persistence");
			vartype = ezxml_attr(var_xml, "type");
			vardims = ezxml_attr(var_xml, "dimensions");
			varunits = ezxml_attr(var_xml, "units");
			vardesc = ezxml_attr(var_xml, "description");
			varname_in_code = ezxml_attr(var_xml, "name_in_code");
			varpackages = ezxml_attr(var_xml, "packages");

			if (varname == NULL) {
				fprintf(stderr,"ERROR: Variable name missing in var_struct %s\n.", structname);
				return 1;
			}

			if(vartype == NULL) {
				fprintf(stderr,"ERROR: Type attribute missing on variable %s in var_struct %s\n.", varname, structname);
				return 1;
			} else if (strcasecmp("logical", vartype) != 0 && strcasecmp("real", vartype) != 0 &&
					strcasecmp("integer", vartype) != 0 && strcasecmp("text", vartype) != 0) {
				fprintf(stderr,"ERROR: Type attribute on variable %s in var_struct %s is not equal to one of logical, real, integer, or text.\n", varname, structname);
				return 1;
			}

			if (vardims == NULL) {
				fprintf(stderr,"ERROR: Dimensions attribute missing for variable %s in var_struct %s.\n", varname, structname);
				return 1;
			} else {
				if (strcasecmp("", vardims) != 0) {
					string = strdup(vardims);
					err_string = check_dimensions(registry, string);
					free(string);

					if(err_string != NULL) {
						fprintf(stderr,"ERROR: Dimension %s on variable %s in var_struct %s not defined.\n", err_string, varname, structname); 
						return 1;
					}
				}
			}

			persistence = PERSISTENT;
			if (varpersistence != NULL) {
				persistence = check_persistence(varpersistence);

				if(persistence == -1){
					fprintf(stderr, "\ton varaible %s in var_struct %s.\n", varname, structname);
					return -1;
				}
			}

			if(varpackages != NULL && persistence == PERSISTENT){
				string = strdup(varpackages);
				err_string = check_packages(registry, string);
				free(string);

				if (err_string != NULL){
					fprintf(stderr, "ERROR: Package %s used on variable %s in var_struct %s is not defined.\n", err_string, varname, structname);
					return 1;
				}
			} else if ( persistence == SCRATCH && varpackages != NULL ) {
				fprintf(stderr, "ERROR: Packages attribute not allowed on scratch variable %s in var_struct %s.\n", varname, structname);
				return -1;
			} else if ( persistence == SCRATCH && varpackages == NULL && structpackages != NULL) {
				fprintf(stderr, "ERROR: Packages attribute inherited from var_struct %s not allowed on scratch var %s in var_struct %s.\n", structname, varname, structname);
				return -1;
			}

		}
	}

	return 0;
}/*}}}*/

int parse_reg_xml(ezxml_t registry, struct namelist **nls, struct dimension ** dims, struct variable ** vars, struct group_list ** groups, struct package ** pkgs, char * modelname, char * corename, char * version)/*{{{*/
{
	struct namelist * nls_ptr, *nls_ptr2;
	struct namelist * nls_chk_ptr;
	struct dimension * dim_ptr, *dim_ptr2;
	struct variable * var_ptr, *var_ptr2;
	struct dimension_list * dimlist_ptr;
	struct dimension * dimlist_cursor;
	struct group_list * grouplist_ptr;
	struct variable_list * vlist_cursor;
	struct package * pkg_ptr;

	ezxml_t dims_xml, dim_xml;
	ezxml_t structs_xml, var_arr_xml, var_xml;
	ezxml_t nmlrecs_xml, nmlopt_xml;
	ezxml_t packages_xml, package_xml;
	ezxml_t streams_xml, stream_xml;

	const char *dimname, *dimunits, *dimdesc, *dimdef;
	const char *nmlrecname, *nmlrecindef;
	const char *nmloptname, *nmlopttype, *nmloptval, *nmloptunits, *nmloptdesc, *nmloptposvals, *nmloptindef;
	const char *structname, *structlevs, *structpackages;
	const char *vararrname, *vararrtype, *vararrdims, *vararrpersistence, *vararrdefaultval, *vararrpackages;
	const char *varname, *varpersistence, *vartype, *vardims, *varunits, *vardesc, *vararrgroup, *varstreams, *vardefaultval, *varpackages;
	const char *packagename, *packagedesc;
	const char *varname_in_code;
	const char *const_model, *const_core, *const_version;
	const char *streamname;

	char dimensions[2048];
	char *dimension_list;
	char dimension_buffer[128];
	char default_value[1024];

	char *string, *tofree, *token;

	NEW_NAMELIST(nls_ptr)
		NEW_DIMENSION(dim_ptr)
		NEW_VARIABLE(var_ptr)
		NEW_GROUP_LIST(grouplist_ptr);
	NEW_PACKAGE(pkg_ptr);
	*nls = nls_ptr;
	*dims = dim_ptr;
	*vars = var_ptr;
	*groups = grouplist_ptr;
	*pkgs = pkg_ptr;

	snprintf(pkg_ptr->name, 1024, "%c", '\0');

	// Get model information
	const_model = ezxml_attr(registry, "model");
	const_core = ezxml_attr(registry, "core");
	const_version = ezxml_attr(registry, "version");

	if(const_model == NULL)
		sprintf(modelname, "MISSING");
	else
		sprintf(modelname, "%s", const_model);

	if(const_core == NULL)
		sprintf(corename, "MISSING");
	else
		sprintf(corename, "%s", const_core);

	if(const_version == NULL)
		sprintf(version, "MISSING");
	else
		sprintf(version, "%s", const_version);

	// Parse Namelist Records
	for (nmlrecs_xml = ezxml_child(registry, "nml_record"); nmlrecs_xml; nmlrecs_xml = nmlrecs_xml->next){
		nmlrecname = ezxml_attr(nmlrecs_xml, "name");
		nmlrecindef = ezxml_attr(nmlrecs_xml, "in_defaults");
		for (nmlopt_xml = ezxml_child(nmlrecs_xml, "nml_option"); nmlopt_xml; nmlopt_xml = nmlopt_xml->next){
			nmloptname = ezxml_attr(nmlopt_xml, "name");
			nmlopttype = ezxml_attr(nmlopt_xml, "type");
			nmloptval = ezxml_attr(nmlopt_xml, "default_value");
			nmloptunits = ezxml_attr(nmlopt_xml, "units");
			nmloptdesc = ezxml_attr(nmlopt_xml, "description");
			nmloptposvals = ezxml_attr(nmlopt_xml, "possible_values");
			nmloptindef = ezxml_attr(nmlopt_xml, "in_defaults");

			snprintf(nls_ptr->record, 1024, "%s", nmlrecname);
			snprintf(nls_ptr->name, 1024, "%s", nmloptname);

			if(strncmp(nmlopttype, "real", 1024) == 0){
				nls_ptr->vtype = REAL;
			} else if(strncmp(nmlopttype, "integer", 1024) == 0){
				nls_ptr->vtype = INTEGER;
			} else if(strncmp(nmlopttype, "logical", 1024) == 0){
				nls_ptr->vtype = LOGICAL;
			} else if(strncmp(nmlopttype, "character", 1024) == 0){
				nls_ptr->vtype = CHARACTER;
			}

			switch(nls_ptr->vtype){
				case REAL:
					nls_ptr->defval.rval = (float)atof(nmloptval);
					break;
				case INTEGER:
					nls_ptr->defval.ival = atoi(nmloptval);
					break;
				case LOGICAL:
					if(strncmp(nmloptval, "true", 1024) ==0){
						nls_ptr->defval.lval = 1;
					} else if (strncmp(nmloptval, "false", 1024) == 0){
						nls_ptr->defval.lval = 0;
					}
					break;
				case CHARACTER:
					snprintf(nls_ptr->defval.cval, 32, "%s", nmloptval);
					break;
			}

			if(nmloptindef == NULL){
				if(nmlrecindef == NULL){
					nls_ptr->write_in_default = 0;
				} else {
					if(strncmp(nmlrecindef, "true", 1024) == 0){
						nls_ptr->write_in_default = 1;
					} else if (strncmp(nmlrecindef, "false", 1024) == 0){
						nls_ptr->write_in_default = 0;
					}
				}
			} else {
				if(strncmp(nmloptindef, "true", 1024) == 0){
					nls_ptr->write_in_default = 1;
				} else if (strncmp(nmloptindef, "false", 1024) == 0){
					nls_ptr->write_in_default = 0;
				}
			}

			NEW_NAMELIST(nls_ptr->next)
				nls_ptr2 = nls_ptr;
			nls_ptr = nls_ptr->next;
		}
	}

	if(nls_ptr2->next) free(nls_ptr2->next);
	nls_ptr2->next = NULL;

	// Parse Packages
	for (packages_xml = ezxml_child(registry, "packages"); packages_xml; packages_xml = packages_xml->next){
		for (package_xml = ezxml_child(packages_xml, "package"); package_xml; package_xml = package_xml->next){
			packagename = ezxml_attr(package_xml, "name");
			packagedesc = ezxml_attr(package_xml, "description");

			if (strlen(pkg_ptr->name) == 0) {
				snprintf(pkg_ptr->name, 1024, "%s", packagename);
			} else {
				NEW_PACKAGE(pkg_ptr->next);
				pkg_ptr = pkg_ptr->next;
				snprintf(pkg_ptr->name, 1024, "%s", packagename);
			}
		}
	}

	// Parse Dimensions
	for (dims_xml = ezxml_child(registry, "dims"); dims_xml; dims_xml = dims_xml->next){
		for (dim_xml = ezxml_child(dims_xml, "dim"); dim_xml; dim_xml = dim_xml->next){
			dimname = ezxml_attr(dim_xml, "name");
			dimdef = ezxml_attr(dim_xml, "definition");	
			dimunits = ezxml_attr(dim_xml, "units");
			dimdesc = ezxml_attr(dim_xml, "description");

			dim_ptr->namelist_defined = 0;

			snprintf(dim_ptr->name_in_file, 1024, "%s", dimname);
			if(dimdef == NULL){
				snprintf(dim_ptr->name_in_code, 1024, "%s", dimname);
				dim_ptr->constant_value = -1;
			} else {
				snprintf(dim_ptr->name_in_code, 1024, "%s", dimdef);
				// Check namelist defined ??
				dim_ptr->constant_value = is_integer_constant(dim_ptr->name_in_code);
				if(strncmp(dim_ptr->name_in_code, "namelist:", 9) == 0) {
					dim_ptr->namelist_defined = 1;
					snprintf(dim_ptr->name_in_code, 1024, "%s", (dim_ptr->name_in_code)+9);

					/* Check that the referenced namelist variable is defined as an integer variable */
					nls_chk_ptr = (*nls)->next;
					while (nls_chk_ptr) {
						if (strncmp(nls_chk_ptr->name, dim_ptr->name_in_code, 1024) == 0) {
							if (nls_chk_ptr->vtype != INTEGER) {
								fprintf(stderr,"\nRegistry error: Namelist variable %s must be an integer for namelist-derived dimension %s\n\n", nls_chk_ptr->name, dim_ptr->name_in_file);
								return 1;
							}
							break;
						} 
						nls_chk_ptr = nls_chk_ptr->next;
					}
					if (!nls_chk_ptr) {
						fprintf(stderr,"\nRegistry error: Namelist variable %s not defined for namelist-derived dimension %s\n\n", dim_ptr->name_in_code, dim_ptr->name_in_file);
						return 1;
					}

				}
			}

			NEW_DIMENSION(dim_ptr->next)
				dim_ptr2 = dim_ptr;
			dim_ptr = dim_ptr->next;
		}   
	}

	if(dim_ptr2->next) free(dim_ptr2->next);
	dim_ptr2->next = NULL;

	// Parse Variable Structures
	for(structs_xml = ezxml_child(registry, "var_struct"); structs_xml; structs_xml = structs_xml->next){
		structname = ezxml_attr(structs_xml, "name");
		structlevs = ezxml_attr(structs_xml, "time_levs");
		structpackages = ezxml_attr(structs_xml, "packages");

		grouplist_ptr = *groups;
		while(grouplist_ptr->next) grouplist_ptr = grouplist_ptr->next;
		NEW_GROUP_LIST(grouplist_ptr->next);
		grouplist_ptr = grouplist_ptr->next;
		snprintf(grouplist_ptr->name, 1024, "%s", structname);
		grouplist_ptr->ntime_levs = atoi(structlevs);
		vlist_cursor = NULL;

		// Parse variable arrays
		for(var_arr_xml = ezxml_child(structs_xml, "var_array"); var_arr_xml; var_arr_xml = var_arr_xml->next){
			vararrname = ezxml_attr(var_arr_xml, "name");
			vararrtype = ezxml_attr(var_arr_xml, "type");
			vararrdims = ezxml_attr(var_arr_xml, "dimensions");
			vararrpersistence = ezxml_attr(var_arr_xml, "persistence");
			vararrdefaultval = ezxml_attr(var_arr_xml, "default_value");
			vararrpackages = ezxml_attr(var_arr_xml, "packages");

			//Parse variables in variable arrays
			for(var_xml = ezxml_child(var_arr_xml, "var"); var_xml; var_xml = var_xml->next){
				varname = ezxml_attr(var_xml, "name");
				varunits = ezxml_attr(var_xml, "units");
				vardesc = ezxml_attr(var_xml, "description");
				vararrgroup = ezxml_attr(var_xml, "array_group");
				varname_in_code = ezxml_attr(var_xml, "name_in_code");
				varpackages = ezxml_attr(var_xml, "packages");

				if(vlist_cursor == NULL){
					NEW_VARIABLE_LIST(grouplist_ptr->vlist);
					vlist_cursor = grouplist_ptr->vlist;
				} else {
					NEW_VARIABLE_LIST(vlist_cursor->next);
					vlist_cursor->next->prev = vlist_cursor;
					vlist_cursor = vlist_cursor->next;
				}
				vlist_cursor->var = var_ptr;
				vlist_cursor->next = NULL;

				var_ptr->ndims = 0;
				var_ptr->timedim = 0;
				var_ptr->iostreams = 0;

				snprintf(var_ptr->name_in_file, 1024, "%s", varname);

				var_ptr->persistence = PERSISTENT;
				if(vararrpersistence != NULL){
					var_ptr->persistence = check_persistence(vararrpersistence);
					if (var_ptr->persistence == -1) return 1;
				}

				/* Check var_arr packages attribute */
				if(varpackages == NULL) {
					varpackages = ezxml_attr(var_arr_xml, "packages");
				}

				/* Check var_struct packages attribute */
				if(varpackages == NULL) {
					varpackages = ezxml_attr(structs_xml, "packages");
				}

				if(varpackages != NULL && var_ptr->persistence == PERSISTENT){
					var_ptr->persistence = PACKAGE;
				}

				if(var_ptr->persistence == PACKAGE) {
					NEW_PACKAGE(var_ptr->package_name);

					string = strdup(varpackages);
					tofree = string;
					token = strsep(&string, ";");

					snprintf(var_ptr->package_name->name, 1024, "%s", token);
					pkg_ptr = var_ptr->package_name;

					while( (token = strsep(&string, ";")) != NULL) {
						NEW_PACKAGE(pkg_ptr->next);						
						pkg_ptr = pkg_ptr->next;
						snprintf(pkg_ptr->name, 1024, "%s", token);
					}
				}

				if(strncmp(vararrtype, "real", 1024) == 0){
					var_ptr->vtype = REAL;
					snprintf(default_value, 1024, "0.0_RKIND");
				} else if(strncmp(vararrtype, "integer", 1024) == 0){
					var_ptr->vtype = INTEGER;
					snprintf(default_value, 1024, "0");
				} else if(strncmp(vararrtype, "logical", 1024) == 0){
					var_ptr->vtype = LOGICAL;
					snprintf(default_value, 1024, ".false.");
				} else if(strncmp(vararrtype, "text", 1024) == 0){
					var_ptr->vtype = CHARACTER;
					snprintf(default_value, 1024, "''");
				}

				NEW_DIMENSION_LIST(dimlist_ptr)
					var_ptr->dimlist = dimlist_ptr;

				snprintf(dimensions,2048, "%s", vararrdims);
				dimension_list = strtok(dimensions, " ");
				while(dimension_list != NULL){
					snprintf(dimension_buffer, 128, "%s", dimension_list);
					if(strncmp(dimension_buffer, "Time", 1024) == 0){
						var_ptr->timedim = 1;
					} else {
						NEW_DIMENSION_LIST(dimlist_ptr->next)
							dimlist_ptr->next->prev = dimlist_ptr;
						dimlist_ptr = dimlist_ptr->next;

						dimlist_cursor = (*dims);
						while(dimlist_cursor && (strncmp(dimension_buffer, dimlist_cursor->name_in_file, 1024) != 0)){
							dimlist_cursor = dimlist_cursor->next;
						}
						if (dimlist_cursor) {
							dimlist_ptr->dim = dimlist_cursor;
						} else {
							fprintf(stderr, "Error: Unknown dimension %s for variable %s\n", dimension_buffer, var_ptr->name_in_file);
							return 1;
						}
						var_ptr->ndims++;
					}
					dimension_list = strtok(NULL, " ");
				}
				dimlist_ptr = var_ptr->dimlist;
				if(var_ptr->dimlist) var_ptr->dimlist = var_ptr->dimlist->next;
				free(dimlist_ptr);

				if(varname_in_code == NULL){
					snprintf(var_ptr->name_in_code, 1024, "%s", varname);
				} else {
					snprintf(var_ptr->name_in_code, 1024, "%s", varname_in_code);
				}

				if(vararrdefaultval == NULL){
					snprintf(var_ptr->default_value, 1024, "%s", default_value);
				} else {
					snprintf(var_ptr->default_value, 1024, "%s", vararrdefaultval);
				}

				snprintf(var_ptr->var_array, 1024, "%s", vararrname);
				snprintf(var_ptr->array_class, 1024, "%s", vararrgroup);

				NEW_VARIABLE(var_ptr->next);
				var_ptr2 = var_ptr;
				var_ptr = var_ptr->next;
			}
		}

		for(var_xml = ezxml_child(structs_xml, "var"); var_xml; var_xml = var_xml->next){
			varname = ezxml_attr(var_xml, "name");
			varpersistence = ezxml_attr(var_xml, "persistence");
			vartype = ezxml_attr(var_xml, "type");
			vardims = ezxml_attr(var_xml, "dimensions");
			varunits = ezxml_attr(var_xml, "units");
			vardesc = ezxml_attr(var_xml, "description");
			varname_in_code = ezxml_attr(var_xml, "name_in_code");
			vardefaultval = ezxml_attr(var_xml, "default_value");
			varpackages = ezxml_attr(var_xml, "packages");

			if(vlist_cursor == NULL){
				NEW_VARIABLE_LIST(grouplist_ptr->vlist);
				vlist_cursor = grouplist_ptr->vlist;
			} else {
				NEW_VARIABLE_LIST(vlist_cursor->next);
				vlist_cursor->next->prev = vlist_cursor;
				vlist_cursor = vlist_cursor->next;
			}
			vlist_cursor->var = var_ptr;
			vlist_cursor->next = NULL;

			var_ptr->ndims = 0;
			var_ptr->timedim = 0;
			var_ptr->iostreams = 0;

			snprintf(var_ptr->name_in_file, 1024, "%s", varname);

			var_ptr->persistence = PERSISTENT;
			if(varpersistence != NULL){
				var_ptr->persistence = check_persistence(varpersistence);
				if(var_ptr->persistence == -1) return 1;
			}

			/* Check packages attribute on var_struct */
			if(varpackages == NULL){
				varpackages = ezxml_attr(structs_xml, "packages");
			}

			if(varpackages != NULL && var_ptr->persistence == PERSISTENT){
				var_ptr->persistence = PACKAGE;
			}

			if(var_ptr->persistence == PACKAGE) {
				NEW_PACKAGE(var_ptr->package_name);

				string = strdup(varpackages);
				tofree = string;
				token = strsep(&string, ";");

				snprintf(var_ptr->package_name->name, 1024, "%s", token);
				pkg_ptr = var_ptr->package_name;

				while( (token = strsep(&string, ";")) != NULL) {
					NEW_PACKAGE(pkg_ptr->next);						
					pkg_ptr = pkg_ptr->next;
					snprintf(pkg_ptr->name, 1024, "%s", token);
				}
			}


			if(strncmp(vartype, "real", 1024) == 0){
				var_ptr->vtype = REAL;
				snprintf(default_value, 1024, "0.0_RKIND");
			} else if(strncmp(vartype, "integer", 1024) == 0){
				var_ptr->vtype = INTEGER;
				snprintf(default_value, 1024, "0");
			} else if(strncmp(vartype, "logical", 1024) == 0){
				var_ptr->vtype = LOGICAL;
				snprintf(default_value, 1024, ".false.");
			} else if(strncmp(vartype, "text", 1024) == 0){
				var_ptr->vtype = CHARACTER;
				snprintf(default_value, 1024, "''");
			}

			NEW_DIMENSION_LIST(dimlist_ptr)
				var_ptr->dimlist = dimlist_ptr;

			snprintf(dimensions, 2048, "%s", vardims);
			dimension_list = strtok(dimensions, " ");
			while(dimension_list != NULL){
				snprintf(dimension_buffer, 128, "%s", dimension_list);
				if(strncmp(dimension_buffer, "Time", 1024) == 0){
					var_ptr->timedim = 1;
				} else {
					NEW_DIMENSION_LIST(dimlist_ptr->next)
						dimlist_ptr->next->prev = dimlist_ptr;
					dimlist_ptr = dimlist_ptr->next;

					dimlist_cursor = (*dims);
					while(dimlist_cursor && (strncmp(dimension_buffer, dimlist_cursor->name_in_file, 1024) != 0) )
						dimlist_cursor = dimlist_cursor->next;
					if (dimlist_cursor) {
						dimlist_ptr->dim = dimlist_cursor;
					} else {
						fprintf(stderr, "Error: Unknown dimension %s for variable %s\n", dimension_buffer, var_ptr->name_in_file);
						return 1;
					}
					var_ptr->ndims++;
				}
				dimension_list = strtok(NULL, " ");
			}

			dimlist_ptr = var_ptr->dimlist;
			if(var_ptr->dimlist) var_ptr->dimlist = var_ptr->dimlist->next;
			free(dimlist_ptr);


			if(varname_in_code == NULL){
				snprintf(var_ptr->name_in_code, 1024, "%s", varname);
			} else {
				snprintf(var_ptr->name_in_code, 1024, "%s", varname_in_code);
			}

			snprintf(var_ptr->var_array, 1024, "-");
			snprintf(var_ptr->array_class, 1024, "-");

			if(vardefaultval == NULL){
				snprintf(var_ptr->default_value, 1024, "%s", default_value);
			} else {
				snprintf(var_ptr->default_value, 1024, "%s", vardefaultval);
			}

			NEW_VARIABLE(var_ptr->next);
			var_ptr2 = var_ptr;
			var_ptr = var_ptr->next;
		}
	}

	if(var_ptr2->next) free(var_ptr2->next);
	var_ptr2->next = NULL;

	grouplist_ptr = *groups;
	if ((*groups)->next) *groups = (*groups)->next;
	if (grouplist_ptr) free(grouplist_ptr);

	// Parse streams
	for (streams_xml = ezxml_child(registry, "streams"); streams_xml; streams_xml = streams_xml->next) {
		for (stream_xml = ezxml_child(streams_xml, "stream"); stream_xml; stream_xml = stream_xml->next) {
			streamname = ezxml_attr(stream_xml, "name");
			if (streamname != NULL) {     /* this should be assured by validate_reg_xml() */
				for (var_xml = ezxml_child(stream_xml, "var"); var_xml; var_xml = var_xml->next) {
					varname = ezxml_attr(var_xml, "name");
					if (varname != NULL) {     /* this should be assured by validate_reg_xml() */
						var_ptr = *vars;
						while (var_ptr) {
							if (strcmp(var_ptr->name_in_file, varname) == 0) {
								if (strcmp(streamname, "output") == 0) {
									var_ptr->iostreams |= OUTPUT0;
								}
								else if (strcmp(streamname, "surface") == 0) {
									var_ptr->iostreams |= SFC0;
								}
								else if (strcmp(streamname, "input") == 0) {
									var_ptr->iostreams |= INPUT0;
								}
								else if (strcmp(streamname, "restart") == 0) {
									var_ptr->iostreams |= RESTART0;
								}
								break;
							}
							var_ptr = var_ptr->next;
						}
					}
					if (var_ptr == NULL) printf("did not find match for %s\n", varname);
				}
			}
		}
	}


	return 0;
}/*}}}*/

int getword(FILE * regfile, char * word)/*{{{*/
{
	int i;
	int c;

	i = 0;

	do { c = getc(regfile); } while (((char)c == ' ' || (char)c == '\n' || (char)c == '\t') && c != EOF);

	while ((char)c == '%') {
		do { c = getc(regfile); } while ((char)c != '\n' && c != EOF);
		do { c = getc(regfile); } while (((char)c == ' ' || (char)c == '\n' || (char)c == '\t') && c != EOF);
	};
	while((char)c != ' ' && (char)c != '\n' && (char)c != '\t' && c != EOF && (char)c != '%') {
		word[i++] = (char)c; 
		c = (char)getc(regfile);
	} 
	word[i] = '\0';

	if ((char)c == '%') do { c = getc(regfile); } while ((char)c != '\n' && c != EOF);

	fprintf(stdout,"%s ",word);
	return c;
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

void sort_vars(struct variable * vars)/*{{{*/
{
	struct variable * var_ptr;
	struct variable * var_ptr2;
	struct variable * var_ptr2_prev;
	char var_array[1024];
	char array_class[1024];

	var_ptr = vars;

	/* Attempt at sorting first on super-array, then on class in the same loop
	   while (var_ptr) {
	   memcpy(var_array, var_ptr->var_array, 1024);
	   memcpy(array_class, var_ptr->array_class, 1024);
	   var_ptr2_prev = var_ptr;
	   var_ptr2 = var_ptr->next;
	   if (var_ptr2 && 
	   (strncmp(var_array, var_ptr2->var_array, 1024) != 0 || strncmp(array_class, var_ptr2->array_class, 1024) != 0)) {
	   while (var_ptr2) {
	   if (strncmp(var_array, var_ptr2->var_array, 1024) == 0 && strncmp(array_class, var_ptr2->array_class, 1024) == 0) {
	   var_ptr2_prev->next = var_ptr2->next;
	   var_ptr2->next = var_ptr->next;
	   var_ptr->next = var_ptr2;
	   var_ptr2 = var_ptr2_prev->next;
	   }
	   else {
	   var_ptr2_prev = var_ptr2_prev->next;
	   var_ptr2 = var_ptr2->next;
	   }
	   }
	   } 
	   var_ptr = var_ptr->next;
	   }
	   */

	while (var_ptr) {
		memcpy(var_array, var_ptr->var_array, 1024);
		var_ptr2_prev = var_ptr;
		var_ptr2 = var_ptr->next;
		if (var_ptr2 && strncmp(var_array, var_ptr2->var_array, 1024) != 0) {
			while (var_ptr2) {
				if (strncmp(var_array, var_ptr2->var_array, 1024) == 0) {
					var_ptr2_prev->next = var_ptr2->next;
					var_ptr2->next = var_ptr->next;
					var_ptr->next = var_ptr2;
					var_ptr2 = var_ptr2_prev->next;
				}
				else {
					var_ptr2_prev = var_ptr2_prev->next;
					var_ptr2 = var_ptr2->next;
				}
			}
		} 
		var_ptr = var_ptr->next;
	}

	var_ptr = vars;

	while (var_ptr) {
		memcpy(array_class, var_ptr->array_class, 1024);
		var_ptr2_prev = var_ptr;
		var_ptr2 = var_ptr->next;
		if (var_ptr2 && strncmp(array_class, var_ptr2->array_class, 1024) != 0) {
			while (var_ptr2) {
				if (strncmp(array_class, var_ptr2->array_class, 1024) == 0) {
					var_ptr2_prev->next = var_ptr2->next;
					var_ptr2->next = var_ptr->next;
					var_ptr->next = var_ptr2;
					var_ptr2 = var_ptr2_prev->next;
				}
				else {
					var_ptr2_prev = var_ptr2_prev->next;
					var_ptr2 = var_ptr2->next;
				}
			}
		} 
		var_ptr = var_ptr->next;
	}
}/*}}}*/

void sort_group_vars(struct group_list * groups)/*{{{*/
{
	struct variable_list * var_list;
	struct variable_list * var_ptr;
	struct variable_list * var_ptr2;
	struct variable_list * var_ptr2_prev;
	struct group_list * group_ptr;
	char var_array[1024];
	char array_class[1024];

	group_ptr = groups;

	while (group_ptr) {

		var_ptr = group_ptr->vlist;

		while (var_ptr) {
			memcpy(var_array, var_ptr->var->var_array, 1024);
			var_ptr2_prev = var_ptr;
			var_ptr2 = var_ptr->next;
			if (var_ptr2 != NULL && strncmp(var_array, var_ptr2->var->var_array, 1024) != 0) {
				while (var_ptr2) {
					if (strncmp(var_array, var_ptr2->var->var_array, 1024) == 0) {
						var_ptr2_prev->next = var_ptr2->next;
						var_ptr2->next = var_ptr->next;
						var_ptr->next = var_ptr2;
						var_ptr2 = var_ptr2_prev->next;
					}
					else {
						var_ptr2_prev = var_ptr2_prev->next;
						var_ptr2 = var_ptr2->next;
					}
				}
			} 
			var_ptr = var_ptr->next;
		}

		var_ptr = group_ptr->vlist;

		while (var_ptr) {
			memcpy(array_class, var_ptr->var->array_class, 1024);
			var_ptr2_prev = var_ptr;
			var_ptr2 = var_ptr->next;
			if (var_ptr2 && strncmp(array_class, var_ptr2->var->array_class, 1024) != 0) {
				while (var_ptr2) {
					if (strncmp(array_class, var_ptr2->var->array_class, 1024) == 0) {
						var_ptr2_prev->next = var_ptr2->next;
						var_ptr2->next = var_ptr->next;
						var_ptr->next = var_ptr2;
						var_ptr2 = var_ptr2_prev->next;
					}
					else {
						var_ptr2_prev = var_ptr2_prev->next;
						var_ptr2 = var_ptr2->next;
					}
				}
			} 
			var_ptr = var_ptr->next;
		}

		group_ptr = group_ptr->next;
	}
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

char * check_streams(char * streams){/*{{{*/
	char * stream;
	int length, i, bad_streams;

	length = strlen(streams);

	stream = (char *)malloc(2*sizeof(char));
	stream[1] = '\0';

	for (i = 0; i < length; i++){
		bad_streams = 1;	
		stream[0] = streams[i];
		if(strcmp(stream, "i") == 0 || strcmp(stream, "r") == 0 || strcmp(stream, "o") == 0 || strcmp(stream, "s") == 0){
			bad_streams = 0;
		}

		if (bad_streams == 1){
			return stream;
		}
	}

	return NULL;
}/*}}}*/

int check_persistence(const char * persistence){/*{{{*/
	if(strncmp(persistence, "persistent", 1024) == 0){
		return PERSISTENT;
	} else if(strncmp(persistence, "scratch", 1024) == 0){
		return SCRATCH;
	} else {
		fprintf(stderr, "ERROR: In check_persistence. Persistence not equal to persistent or scratch.\n");
		return -1;
	}
}/*}}}*/
