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
#include "ezxml.h"
#include "fortprintf.h"
#include "registry_types.h"
#include "gen_inc.h"
#include "utility.h"


int is_unique_field(ezxml_t registry, ezxml_t field, const char *check_name);
int is_unique_struct(ezxml_t registry, ezxml_t check_struct, const char *check_name);
int check_for_unique_names(ezxml_t registry, ezxml_t current_position);
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

	return 0;
}/*}}}*/


int validate_reg_xml(ezxml_t registry)/*{{{*/
{
	ezxml_t dims_xml, dim_xml;
	ezxml_t structs_xml, var_arr_xml, var_xml, stream_var_xml;
	ezxml_t nmlrecs_xml, nmlopt_xml;
	ezxml_t streams_xml, stream_xml, substream_xml;
	ezxml_t streams_xml2, stream_xml2;

	const char *dimname, *dimunits, *dimdesc, *dimdef, *dimdecomp;
	const char *nmlrecname, *nmlrecindef;
	const char *nmloptname, *nmlopttype, *nmloptval, *nmloptunits, *nmloptdesc, *nmloptposvals, *nmloptindef;
	const char *nmloptbounds, *nmloptcellmeasures, *nmloptcellmethod, *nmloptcoords, *nmloptstdname;
	const char *structname, *structpackages, *structstreams;
	const char *vararrname, *vararrtype, *vararrdims, *vararrpersistence, *vararrpackages, *vararrstreams;
	const char *vararrbounds, *vararrcellmeasures, *vararrcellmethod, *vararrcoords, *vararrstdname;
	const char *varname, *varpersistence, *vartype, *vardims, *varunits, *vardesc, *vararrgroup, *varstreams, *varpackages;
	const char *varbounds, *varcellmeasures, *varcellmethod, *varcoords, *varstdname;
	const char *varname_in_code, *varname_in_stream;
	const char *const_model, *const_core, *const_version;
	const char *streamname, *streamtype, *streamfilename, *streamrecords, *streaminterval_in, *streaminterval_out, *streampackages, *streamvarpackages;
	const char *streamimmutable, *streamformat;
	const char *substreamname, *streamimmutable2;
	const char *streamname2, *streamtype2, *streamfilename2;
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
			nmloptbounds = ezxml_attr(nmlopt_xml, "bounds");
			nmloptcellmeasures = ezxml_attr(nmlopt_xml, "cell_measures");
			nmloptcellmethod = ezxml_attr(nmlopt_xml, "cell_method");
			nmloptcoords = ezxml_attr(nmlopt_xml, "coordinates");
			nmloptstdname = ezxml_attr(nmlopt_xml, "standard_name");
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
			dimdecomp = ezxml_attr(dim_xml, "decomposition");

			if (dimname == NULL){
				fprintf(stderr,"ERROR: Name missing for dimension.\n");
				return 1;
			}

			if (dimdef != NULL){
				if ( dimdecomp != NULL ) {
					fprintf(stderr, "ERROR: Dimension %s cannot have a decomposition and a definition attribute.\n", dimname);
					return 1;
				}
				if (strncmp(dimdef, "namelist:", 9) == 0){
					found = 0;
					snprintf(name_holder, 1024, "%s",dimdef);
					snprintf(name_holder, 1024, "%s",(name_holder)+9);
					for (nmlrecs_xml = ezxml_child(registry, "nml_record"); nmlrecs_xml; nmlrecs_xml = nmlrecs_xml->next){

						for (nmlopt_xml = ezxml_child(nmlrecs_xml, "nml_option"); nmlopt_xml; nmlopt_xml = nmlopt_xml->next){
							nmloptname = ezxml_attr(nmlopt_xml, "name");
							nmlopttype = ezxml_attr(nmlopt_xml, "type");

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
		structstreams = ezxml_attr(structs_xml, "streams");

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
				fprintf(stderr, "ERROR: Package %s attached to var_struct %s is not defined.\n", err_string, structname);
				return 1;
			}
		}

		if (structstreams != NULL) {
			string = strdup(structstreams);
			err_string = check_streams(registry, string);
			free(string);

			if (err_string != NULL) {
				fprintf(stderr, "ERROR: Stream %s attached to var_struct %s is not defined.\n", err_string, structname);
				return 1;
			}
		}

		// Validate variable arrays
		for(var_arr_xml = ezxml_child(structs_xml, "var_array"); var_arr_xml; var_arr_xml = var_arr_xml->next){
			vararrname = ezxml_attr(var_arr_xml, "name");
			vararrtype = ezxml_attr(var_arr_xml, "type");
			vararrdims = ezxml_attr(var_arr_xml, "dimensions");
			vararrpersistence = ezxml_attr(var_arr_xml, "persistence");
                        vararrbounds = ezxml_attr(var_arr_xml, "bounds");
                        vararrcellmeasures = ezxml_attr(var_arr_xml, "cell_measures");
                        vararrcellmethod = ezxml_attr(var_arr_xml, "cell_method");
                        vararrcoords = ezxml_attr(var_arr_xml, "coordinates");
                        vararrstdname = ezxml_attr(var_arr_xml, "standard_name");
			vararrpackages = ezxml_attr(var_arr_xml, "packages");
			vararrstreams = ezxml_attr(var_arr_xml, "streams");
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
					fprintf(stderr, "ERROR: Package %s attached to var_array %s in var_struct %s is not defined.\n", err_string, vararrname, structname);
					return 1;
				}
			}

			if (persistence == SCRATCH && vararrstreams != NULL){
				fprintf(stderr, "ERROR: Streams attribute not allowed on scratch var_array %s in var_struct %s.\n", vararrname, structname);
				return -1;
			} 
			else if (persistence == SCRATCH && vararrstreams == NULL && structstreams != NULL) {
				fprintf(stderr, "ERROR: Streams attribute inherited from var_struct %s not allowed on scratch var_array %s in var_struct %s.\n", structname, vararrname, structname);
				return -1;
			} 
			else if (persistence == PERSISTENT && vararrstreams != NULL) {
				string = strdup(vararrstreams);
				err_string = check_streams(registry, string);
				free(string);

				if (err_string != NULL) {
					fprintf(stderr, "ERROR: Stream %s attached to var_array %s in var_struct %s is not defined.\n", err_string, vararrname, structname);
					return 1;
				}
			}


			// Validate variables in variable arrays
			for(var_xml = ezxml_child(var_arr_xml, "var"); var_xml; var_xml = var_xml->next){
				varname = ezxml_attr(var_xml, "name");
				varunits = ezxml_attr(var_xml, "units");
				varbounds = ezxml_attr(var_xml, "bounds");
				varcellmeasures = ezxml_attr(var_xml, "cell_measures");
				varcellmethod = ezxml_attr(var_xml, "cell_method");
				varcoords = ezxml_attr(var_xml, "coordinates");
				varstdname = ezxml_attr(var_xml, "standard_name");
				vardesc = ezxml_attr(var_xml, "description");
				vararrgroup = ezxml_attr(var_xml, "array_group");
				varname_in_code = ezxml_attr(var_xml, "name_in_code");
				varpackages = ezxml_attr(var_xml, "packages");
				varstreams = ezxml_attr(var_xml, "streams");

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

				if (persistence == SCRATCH && vararrstreams != NULL) {
					fprintf(stderr, "ERROR: Streams attribute not allowed on constituent variable %s within scratch var_srray %s in var_struct %s.\n", varname, vararrname, structname);
					return 1;
				}

				if(varpackages != NULL){
					string = strdup(varpackages);
					err_string = check_packages(registry, string);
					free(string);

					if (err_string != NULL){
						fprintf(stderr, "ERROR: Package %s attached to constituent variable %s in var_array %s var_struct %s is not defined.\n", err_string, varname, vararrname, structname);
						return 1;
					}
				}

				if(varstreams != NULL){
					string = strdup(varstreams);
					err_string = check_streams(registry, string);
					free(string);

					if (err_string != NULL){
						fprintf(stderr, "ERROR: Stream %s attached to constituent variable %s in var_array %s var_struct %s is not defined.\n", err_string, varname, vararrname, structname);
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
			varbounds = ezxml_attr(var_xml, "bounds");
			varcellmeasures = ezxml_attr(var_xml, "cell_measures");
			varcellmethod = ezxml_attr(var_xml, "cell_method");
			varcoords = ezxml_attr(var_xml, "coordinates");
			varstdname = ezxml_attr(var_xml, "standard_name");
			vardesc = ezxml_attr(var_xml, "description");
			varname_in_code = ezxml_attr(var_xml, "name_in_code");
			varpackages = ezxml_attr(var_xml, "packages");
			varstreams = ezxml_attr(var_xml, "streams");
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
					fprintf(stderr, "ERROR: Package %s attached to variable %s in var_struct %s is not defined.\n", err_string, varname, structname);
					return 1;
				}
			} else if ( persistence == SCRATCH && varpackages != NULL ) {
				fprintf(stderr, "ERROR: Packages attribute not allowed on scratch variable %s in var_struct %s.\n", varname, structname);
				return -1;
			} else if ( persistence == SCRATCH && varpackages == NULL && structpackages != NULL) {
				fprintf(stderr, "ERROR: Packages attribute inherited from var_struct %s not allowed on scratch var %s in var_struct %s.\n", structname, varname, structname);
				return -1;
			}

			if (varstreams != NULL && persistence == PERSISTENT) {
				string = strdup(varstreams);
				err_string = check_streams(registry, string);
				free(string);

				if (err_string != NULL) {
					fprintf(stderr, "ERROR: Stream %s attached to variable %s in var_struct %s is not defined.\n", err_string, varname, structname);
					return 1;
				}
			} 
			else if ( persistence == SCRATCH && varstreams != NULL ) {
				fprintf(stderr, "ERROR: Streams attribute not allowed on scratch variable %s in var_struct %s.\n", varname, structname);
				return -1;
			} 
			else if ( persistence == SCRATCH && varstreams == NULL && structstreams != NULL) {
				fprintf(stderr, "ERROR: Streams attribute inherited from var_struct %s not allowed on scratch var %s in var_struct %s.\n", structname, varname, structname);
				return -1;
			}

		}
	}

	// Validate default streams
	for (streams_xml = ezxml_child(registry, "streams"); streams_xml; streams_xml = streams_xml->next) {
		for (stream_xml = ezxml_child(streams_xml, "stream"); stream_xml; stream_xml = stream_xml->next) {
			streamname = ezxml_attr(stream_xml, "name");
			streamtype = ezxml_attr(stream_xml, "type");
			streamfilename = ezxml_attr(stream_xml, "filename_template");
			streaminterval_in = ezxml_attr(stream_xml, "input_interval");
			streaminterval_out = ezxml_attr(stream_xml, "output_interval");
			streampackages = ezxml_attr(stream_xml, "packages");
			streamimmutable = ezxml_attr(stream_xml, "immutable");
			streamformat = ezxml_attr(stream_xml, "runtime_format");

			if (streamname == NULL) {
				fprintf(stderr, "ERROR: Stream specification missing \"name\" attribute.\n");
				return 1;
			}
			else if (streamtype == NULL) {
				fprintf(stderr, "ERROR: Stream specification for %s missing \"type\" attribute.\n", streamname);
				return 1;
			}
			else if (streamfilename == NULL) {
				fprintf(stderr, "ERROR: Stream specification for %s missing \"filename_template\" attribute.\n", streamname);
				return 1;
			}
			else if (strstr(streamtype, "input") != NULL && streaminterval_in == NULL) {
				fprintf(stderr, "ERROR: Stream %s is marked as input but is missing \"input_interval\" attribute.\n", streamname);
				return 1;
			}
			else if (strstr(streamtype, "output") != NULL && streaminterval_out == NULL) {
				fprintf(stderr, "ERROR: Stream %s is marked as output but is missing \"output_interval\" attribute.\n", streamname);
				return 1;
			}
			else if (streamformat == NULL && (streamimmutable == NULL || (streamimmutable != NULL && strcmp(streamimmutable,"true") != 0))) {
				fprintf(stderr, "ERROR: Mutable stream %s must have the \"runtime_format\" attribute.\n", streamname);
				return 1;
			}
			else {
				/* Check that each stream added to an immutable stream is immutable */
				if (streamimmutable != NULL && strcmp(streamimmutable, "true") == 0) {
					for (substream_xml = ezxml_child(stream_xml, "stream"); substream_xml; substream_xml = substream_xml->next){
						substreamname = ezxml_attr(substream_xml, "name");
						found = 0;

						for (streams_xml2 = ezxml_child(registry, "streams"); streams_xml2; streams_xml2 = streams_xml2->next){
							for (stream_xml2 = ezxml_child(streams_xml2, "stream"); stream_xml2; stream_xml2 = stream_xml2->next){
								streamname2 = ezxml_attr(stream_xml2, "name");

								if (substreamname != NULL && streamname2 != NULL && strcmp(substreamname, streamname2) == 0){
									streamimmutable2 = ezxml_attr(stream_xml2, "immutable");
									found = 1;

									if (streamimmutable2 == NULL || strcmp(streamimmutable2, "true") != 0){
										fprintf(stderr, "ERROR: Immutable stream %s cannot contain mutable streams (e.g. %s).\n", streamname, substreamname);
										return 1;
									}
								}
							}
						}
					}
				}
				for (stream_var_xml = ezxml_child(stream_xml, "var"); stream_var_xml; stream_var_xml = stream_var_xml->next) {
					varname_in_stream = ezxml_attr(stream_var_xml, "name");
					streamvarpackages = ezxml_attr(stream_var_xml, "packages");

					if (varname_in_stream == NULL) {
						fprintf(stderr, "ERROR: Variable field in stream \"%s\" specification missing \"name\" attribute.\n", streamname);
						return 1;
					}

					/* Check that runtime_format is a valid option for mutable streams */
					if (streamimmutable == NULL || (streamimmutable != NULL && strcmp(streamimmutable,"true") != 0)) {
						if (strcmp(streamformat, "single_file") != 0 && strcmp(streamformat, "separate_file") != 0) {
							fprintf(stderr, "ERROR: Runtime_format specification for stream \"%s\" must be either \"single_file\" or \"separate_file\".\n", streamname);
							return 1;
						}
					}


					/* Check that the variable being added to the stream has been defined */
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

					/* did we find what we were looking for? */
					if (var_xml == NULL) {
						fprintf(stderr, "ERROR: Trying to add undefined variable %s to stream %s.\n", varname_in_stream, streamname);
						return 1;
					}	

					if (streamvarpackages != NULL) {
						string = strdup(streamvarpackages);
						err_string = check_packages(registry, string);
						free(string);

						if (err_string != NULL) {
							fprintf(stderr, "ERROR: Package \"%s\" attached to var \"%s\" in stream \"%s\" is not defined.\n", err_string, varname_in_stream, streamname);
							return 1;
						}
					}
				}
				
				/* Validate packages for var_struct members of the stream */
				for (stream_var_xml = ezxml_child(stream_xml, "var_struct"); stream_var_xml; stream_var_xml = stream_var_xml->next) {
					varname_in_stream = ezxml_attr(stream_var_xml, "name");
					streamvarpackages = ezxml_attr(stream_var_xml, "packages");

					if (varname_in_stream == NULL) {
						fprintf(stderr, "ERROR: Variable structure in stream \"%s\" specification missing \"name\" attribute.\n", streamname);
						return 1;
					}

					if (streamvarpackages != NULL) {
						string = strdup(streamvarpackages);
						err_string = check_packages(registry, string);
						free(string);

						if (err_string != NULL) {
							fprintf(stderr, "ERROR: Package \"%s\" attached to var_struct \"%s\" in stream \"%s\" is not defined.\n", err_string, varname_in_stream, streamname);
							return 1;
						}
					}
				}
				
				/* Validate packages for var_array members of the stream */
				for (stream_var_xml = ezxml_child(stream_xml, "var_array"); stream_var_xml; stream_var_xml = stream_var_xml->next) {
					varname_in_stream = ezxml_attr(stream_var_xml, "name");
					streamvarpackages = ezxml_attr(stream_var_xml, "packages");

					if (varname_in_stream == NULL) {
						fprintf(stderr, "ERROR: Variable array in stream \"%s\" specification missing \"name\" attribute.\n", streamname);
						return 1;
					}

					if (streamvarpackages != NULL) {
						string = strdup(streamvarpackages);
						err_string = check_packages(registry, string);
						free(string);

						if (err_string != NULL) {
							fprintf(stderr, "ERROR: Package \"%s\" attached to var_array \"%s\" in stream \"%s\" is not defined.\n", err_string, varname_in_stream, streamname);
							return 1;
						}
					}
				}
				
				/* Validate packages for stream members of the stream */
				for (stream_var_xml = ezxml_child(stream_xml, "stream"); stream_var_xml; stream_var_xml = stream_var_xml->next) {
					varname_in_stream = ezxml_attr(stream_var_xml, "name");
					streamvarpackages = ezxml_attr(stream_var_xml, "packages");

					if (varname_in_stream == NULL) {
						fprintf(stderr, "ERROR: Variable array in stream \"%s\" specification missing \"name\" attribute.\n", streamname);
						return 1;
					}

					if (streamvarpackages != NULL) {
						string = strdup(streamvarpackages);
						err_string = check_packages(registry, string);
						free(string);

						if (err_string != NULL) {
							fprintf(stderr, "ERROR: Package \"%s\" attached to stream \"%s\" in stream \"%s\" is not defined.\n", err_string, varname_in_stream, streamname);
							return 1;
						}
					}
				}
			}

			if (streamformat != NULL && streamimmutable != NULL && strcmp(streamimmutable,"true") == 0) {
				fprintf(stderr, "Warning: runtime_format attribute has no effect for immutable stream \"%s\".\n", streamname);
			}

			if (streampackages != NULL) {
				string = strdup(streampackages);
				err_string = check_packages(registry, string);
				free(string);

				if (err_string != NULL){
					fprintf(stderr, "ERROR: Package \"%s\" attached to stream \"%s\" is not defined.\n", err_string, streamname);
					return 1;
				}
			}

		}
	}
	for (streams_xml = ezxml_child(registry, "streams"); streams_xml; streams_xml = streams_xml->next) {
		for (stream_xml = ezxml_child(streams_xml, "stream"); stream_xml; stream_xml = stream_xml->next) {
			streamname = ezxml_attr(stream_xml, "name");
			streamfilename = ezxml_attr(stream_xml, "filename_template");
			streamtype = ezxml_attr(stream_xml, "type");
			
			/* Check that this stream's filename template is unique among all streams */
			for (streams_xml2 = ezxml_child(registry, "streams"); streams_xml2; streams_xml2 = streams_xml2->next) {
				for (stream_xml2 = ezxml_child(streams_xml2, "stream"); stream_xml2; stream_xml2 = stream_xml2->next) {
					streamname2 = ezxml_attr(stream_xml2, "name");
					streamfilename2 = ezxml_attr(stream_xml2, "filename_template");
					streamtype2 = ezxml_attr(stream_xml, "type");

					if (stream_xml != stream_xml2) {
						if (strcmp(streamfilename, streamfilename2) == 0) {
							if ( strstr(streamtype, "output") != NULL || strstr(streamtype2, "output") != NULL ) {
								fprintf(stderr, "ERROR: Streams %s and %s have a conflicting filename template of %s and one or more has a type that contains output.\n", streamname, streamname2, streamfilename);
								return 1;
							}
						}
					}
				}
			}
		}
	}

	if(check_for_unique_names(registry, registry)){
		fprintf(stderr, "ERROR: Structures and Fields are required to have unique names for I/O reasons.\n");
		fprintf(stderr, "       Please fix duplicates in the Registry.xml file.\n");
		fprintf(stderr, "       You may use the name_in_code attribute to give them the same name inside the model,\n");
		fprintf(stderr, "       but the name attribute is required to be unique.\n");
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

	// Generate code to read and write fields
	err = generate_immutable_streams(registry);

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


int is_unique_struct(ezxml_t current_position, ezxml_t check_struct, const char *check_name){/*{{{*/
	ezxml_t struct_xml;

	const char *name;

	int test;


	test = 1;

	for(struct_xml = ezxml_child(current_position, "var_struct"); struct_xml; struct_xml = struct_xml->next){
		name = ezxml_attr(struct_xml, "name");

		if(strcmp(name, check_name) == 0 && struct_xml != check_struct){
			return 0;
		} else {
			test = is_unique_struct(struct_xml, check_struct, check_name);
			if ( !test ) {
				return 0;
			}
		}
	}

	return 1;
}/*}}}*/


int check_for_unique_names(ezxml_t registry, ezxml_t current_position){/*{{{*/
	ezxml_t struct_xml, var_arr_xml, var_xml;

	const char *name;

	for(struct_xml = ezxml_child(current_position, "var_struct"); struct_xml; struct_xml = struct_xml->next){
		name = ezxml_attr(struct_xml, "name");

		if(!is_unique_struct(registry, struct_xml, name)){
			fprintf(stderr, "ERROR: Struct %s is not uniqe.\n", name);
			return 1;
		}

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

		check_for_unique_names(registry, struct_xml);
	}

	return 0;
}/*}}}*/
