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
#include "fortprintf.h"
#include "registry_types.h"
#include "gen_inc.h"
#include "ezxml/ezxml.h"
#include "utility.h"


int is_unique_field(ezxml_t registry, ezxml_t field, const char *check_name);
int check_for_unique_names(ezxml_t registry);
int is_integer_constant(char *);
int parse_reg_xml(ezxml_t registry);
int validate_reg_xml(ezxml_t registry);


int main(int argc, char ** argv)/*{{{*/
{
	FILE * regfile;
	struct namelist * nls;
	struct dimension * dims;
	struct variable * vars;
	struct group_list * groups;
	struct package * pkgs;
	int err;

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

	// Cleanup registry structures
	err = push_attributes(registry);
	err = merge_structs_and_var_arrays(registry);
	err = merge_streams(registry);

	if (validate_reg_xml(registry)) {
		fprintf(stderr, "Validation failed.....\n");
		return 1;
	}

	write_model_variables(registry);

	if (parse_reg_xml(registry)) {
		fprintf(stderr, "Parsing failed.....\n");
		return 1;
	}

	//write_default_namelist(nls, corename);
	write_default_namelist(registry);

	return 0;
}/*}}}*/


int validate_reg_xml(ezxml_t registry)/*{{{*/
{
	ezxml_t dims_xml, dim_xml;
	ezxml_t structs_xml, var_arr_xml, var_xml, stream_var_xml;
	ezxml_t nmlrecs_xml, nmlopt_xml;
	ezxml_t streams_xml, stream_xml;

	const char *dimname, *dimunits, *dimdesc, *dimdef;
	const char *nmlrecname, *nmlrecindef;
	const char *nmloptname, *nmlopttype, *nmloptval, *nmloptunits, *nmloptdesc, *nmloptposvals, *nmloptindef;
	const char *structname, *structpackages;
	const char *vararrname, *vararrtype, *vararrdims, *vararrpersistence, *vararrpackages;
	const char *varname, *varpersistence, *vartype, *vardims, *varunits, *vardesc, *vararrgroup, *varstreams, *varpackages;
	const char *varname_in_code, *varname_in_stream;
	const char *const_model, *const_core, *const_version;
	const char *streamname;
	const char *streamtype;
	const char *time_levs;

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
		time_levs = ezxml_attr(structs_xml, "time_levs");
		structpackages = ezxml_attr(structs_xml, "packages");

		if (structname == NULL){
			fprintf(stderr,"ERROR: Name missing for var_struct.\n");
			return 1;
		}

		if (time_levs == NULL){
			fprintf(stderr,"ERROR: time_levs attribute missing for var_struct %s.\n", structname);
			return 1;
		} else {
			if (atoi(time_levs) == 0){
				fprintf(stderr, "WARNING: time_levs attribute on var_struct %s is 0. It will be replaced with 1.\n", structname);
			} else if (atoi(time_levs) < 1){
				fprintf(stderr, "ERROR: time_levs attribute on var_struct %s is negative.\n", structname);
				return 1;
			}
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
			time_levs = ezxml_attr(var_arr_xml, "time_levs");

			if (vararrname == NULL){
				fprintf(stderr,"ERROR: Name attribute missing for var_array in var_struct %s.\n", structname);
				return 1;
			}

			if (time_levs != NULL){
				if (atoi(time_levs) == 0){
					fprintf(stderr, "WARNING: time_levs attribute on var_array %s in var_struct %s is 0. It will be replaced with 1.\n", vararrname, structname);
				} else if (atoi(time_levs) < 1){
					fprintf(stderr, "ERROR: time_levs attribute on var_array %s in var_struct %s is negative.\n", vararrname, structname);
					return 1;
				}
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
			time_levs = ezxml_attr(var_xml, "time_levs");

			if (varname == NULL) {
				fprintf(stderr,"ERROR: Variable name missing in var_struct %s\n.", structname);
				return 1;
			}

			if (time_levs != NULL){
				if (atoi(time_levs) == 0){
					fprintf(stderr, "WARNING: time_levs attribute on var %s in var_struct %s is 0. It will be replaced with 1.\n", varname, structname);
				} else if (atoi(time_levs) < 1){
					fprintf(stderr, "ERROR: time_levs attribute on var %s in var_struct %s is negative.\n", varname, structname);
					return 1;
				}
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

	// Validate default streams
	for (streams_xml = ezxml_child(registry, "streams"); streams_xml; streams_xml = streams_xml->next) {
		for (stream_xml = ezxml_child(streams_xml, "stream"); stream_xml; stream_xml = stream_xml->next) {
			streamname = ezxml_attr(stream_xml, "name");
			streamtype = ezxml_attr(stream_xml, "type");
			if (streamname == NULL) {
				fprintf(stderr, "ERROR: Stream specification missing \"name\" attribute.\n");
				return 1;
			}
			else if (streamtype == NULL) {
				fprintf(stderr, "ERROR: Stream specification missing \"type\" attribute.\n");
				return 1;
			}
			else {
				for (stream_var_xml = ezxml_child(stream_xml, "var"); stream_var_xml; stream_var_xml = stream_var_xml->next) {
					varname_in_stream = ezxml_attr(stream_var_xml, "name");
					if (varname_in_stream == NULL) {
						fprintf(stderr, "ERROR: Variable field in stream \"%s\" specification missing \"name\" attribute.\n", streamname);
						return 1;
					}


					// Check that the variable being added to the stream has been defined
					for (structs_xml = ezxml_child(registry, "var_struct"); structs_xml; structs_xml = structs_xml->next) {
						for (var_arr_xml = ezxml_child(structs_xml, "var_array"); var_arr_xml; var_arr_xml = var_arr_xml->next) {
							for (var_xml = ezxml_child(var_arr_xml, "var"); var_xml; var_xml = var_xml->next) {
								varname = ezxml_attr(var_xml, "name");
								if (strcmp(varname, varname_in_stream) == 0) {
									goto done_searching;	
								}
							}
						}
						for (var_xml = ezxml_child(structs_xml, "var"); var_xml; var_xml = var_xml->next) {
							varname = ezxml_attr(var_xml, "name");
							if (strcmp(varname, varname_in_stream) == 0) {
								goto done_searching;	
							}
						}
					}

done_searching:

					// did we find what we were looking for?
					if (var_xml == NULL) {
						fprintf(stderr, "ERROR: Trying to add undefined variable %s to stream %s.\n", varname_in_stream, streamname);
						return 1;
					}	


				}
			}
		}
	}

	if(check_for_unique_names(registry)){
		fprintf(stderr, "ERROR: Fields are required to have unique names for I/O reasons.\n");
		fprintf(stderr, "       Please fix duplicates in the Registry.xml file.\n");
		return 1;
	}

	return 0;
}/*}}}*/


int parse_reg_xml(ezxml_t registry)/*{{{*/
{
	ezxml_t dims_xml, dim_xml;
	ezxml_t structs_xml, var_arr_xml, var_xml;
	ezxml_t nmlrecs_xml, nmlopt_xml;
	ezxml_t packages_xml, package_xml;
	ezxml_t streams_xml, stream_xml;

	int err;


	// Parse Packages
	err = parse_packages_from_registry(registry);

	// Parse namelist records
	err = parse_namelist_records_from_registry(registry);

	// Parse dimensions
	err = parse_dimensions_from_registry(registry);

	// Parse variable structures
	err = parse_structs_from_registry(registry);

	// Generate routines to link fields for multiple blocks
	err = generate_field_links(registry);

	// Generate halo exchange and copy routine
	err = generate_field_halo_exchanges_and_copies(registry);

	// Generate halo exchange and copy routine
	err = generate_field_reads_and_writes(registry);

	return 0;
}/*}}}*/


int is_unique_field(ezxml_t registry, ezxml_t field, const char *check_name){/*{{{*/
	ezxml_t struct_xml, var_arr_xml, var_xml;

	const char *name;

	for(struct_xml = ezxml_child(registry, "var_struct"); struct_xml; struct_xml = struct_xml->next){
		for(var_arr_xml = ezxml_child(struct_xml, "var_array"); var_arr_xml; var_arr_xml = var_arr_xml->next){
			for(var_xml = ezxml_child(var_arr_xml, "var"); var_xml; var_xml = var_xml->next){
				name = ezxml_attr(var_xml, "name");

				if(strcmp(name, check_name) == 0 && var_xml != field){
					return 0;
				}
			}
		}
		for(var_xml = ezxml_child(struct_xml, "var"); var_xml; var_xml = var_xml->next){
			name = ezxml_attr(var_xml, "name");

			if(strcmp(name, check_name) == 0 && var_xml != field){
				return 0;
			}
		}
	}

	return 1;
}/*}}}*/


int check_for_unique_names(ezxml_t registry){/*{{{*/
	ezxml_t struct_xml, var_arr_xml, var_xml;

	const char *name;

	for(struct_xml = ezxml_child(registry, "var_struct"); struct_xml; struct_xml = struct_xml->next){
		for(var_arr_xml = ezxml_child(struct_xml, "var_array"); var_arr_xml; var_arr_xml = var_arr_xml->next){
			for(var_xml = ezxml_child(var_arr_xml, "var"); var_xml; var_xml = var_xml->next){
				name = ezxml_attr(var_xml, "name");
				if(!is_unique_field(registry, var_xml, name)){
					fprintf(stderr, "ERROR: Field %s is not unique.\n", name);
					return 1;
				}
			}
		}

		for(var_xml = ezxml_child(struct_xml, "var"); var_xml; var_xml = var_xml->next){
			name = ezxml_attr(var_xml, "name");
			if(!is_unique_field(registry, var_xml, name)){
				fprintf(stderr, "ERROR: Field %s is not unique.\n", name);
				return 1;
			}
		}
	}

	return 0;
}/*}}}*/
