// Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
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
int parse_reg_xml(FILE * regfile, struct namelist **nls, struct dimension ** dims, struct variable ** vars, struct group_list ** groups, char * modelname, char * corename, char * version);

int main(int argc, char ** argv)
{
   FILE * regfile;
   struct namelist * nls;
   struct dimension * dims;
   struct variable * vars;
   struct group_list * groups;

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
  
   if (parse_reg_xml(regfile, &nls, &dims, &vars, &groups, modelname, corename, version)) {
      return 1;
   }
  
   sort_vars(vars);
   sort_group_vars(groups);

   gen_history_attributes(modelname, corename, version);
   gen_namelists(nls);
   gen_field_defs(groups, vars, dims);
   gen_reads(groups, vars, dims);
   gen_writes(groups, vars, dims, nls);

   return 0;
}

int parse_reg_xml(FILE * regfile, struct namelist **nls, struct dimension ** dims, struct variable ** vars, struct group_list ** groups, char * modelname, char * corename, char * version)
{
	struct namelist * nls_ptr, *nls_ptr2;
	struct namelist * nls_chk_ptr;
	struct dimension * dim_ptr, *dim_ptr2;
	struct variable * var_ptr, *var_ptr2;
	struct dimension_list * dimlist_ptr;
	struct dimension * dimlist_cursor;
	struct group_list * grouplist_ptr;
	struct variable_list * vlist_cursor;

	ezxml_t registry = ezxml_parse_fp(regfile);
	ezxml_t dims_xml, dim_xml;
	ezxml_t structs_xml, var_arr_xml, var_xml;
	ezxml_t nmlrecs_xml, nmlopt_xml;

	const char *dimname, *dimunits, *dimdesc, *dimdef;
	const char *nmlrecname, *nmloptname, *nmlopttype, *nmloptval, *nmloptunits, *nmloptdesc, *nmloptposvals;
	const char *structname, *structlevs;
	const char *vararrname, *vararrtype, *vararrdims, *vararrpersistence;
	const char *varname, *varpersistence, *vartype, *vardims, *varunits, *vardesc, *vararrgroup, *varstreams;
	const char *varname_in_code;
	const char *const_model, *const_core, *const_version;

	char dimensions[2048];
	char *dimension_list;
	char dimension_buffer[128];
	char streams_buffer[128];

	NEW_NAMELIST(nls_ptr)
	NEW_DIMENSION(dim_ptr)
	NEW_VARIABLE(var_ptr)
	NEW_GROUP_LIST(grouplist_ptr);
	*nls = nls_ptr;
	*dims = dim_ptr;
	*vars = var_ptr;
	*groups = grouplist_ptr;

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
		for (nmlopt_xml = ezxml_child(nmlrecs_xml, "nml_option"); nmlopt_xml; nmlopt_xml = nmlopt_xml->next){
			nmloptname = ezxml_attr(nmlopt_xml, "name");
			nmlopttype = ezxml_attr(nmlopt_xml, "type");
			nmloptval = ezxml_attr(nmlopt_xml, "default_value");
			nmloptunits = ezxml_attr(nmlopt_xml, "units");
			nmloptdesc = ezxml_attr(nmlopt_xml, "description");
			nmloptposvals = ezxml_attr(nmlopt_xml, "possible_values");

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

			NEW_NAMELIST(nls_ptr->next)
			nls_ptr2 = nls_ptr;
			nls_ptr = nls_ptr->next;
		}
	}

	if(nls_ptr2->next) free(nls_ptr2->next);
	nls_ptr2->next = NULL;

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
								printf("\nRegistry error: Namelist variable %s must be an integer for namelist-derived dimension %s\n\n", nls_chk_ptr->name, dim_ptr->name_in_file);
								return 1;
							}
							break;
						} 
						nls_chk_ptr = nls_chk_ptr->next;
					}
					if (!nls_chk_ptr) {
						printf("\nRegistry error: Namelist variable %s not defined for namelist-derived dimension %s\n\n", dim_ptr->name_in_code, dim_ptr->name_in_file);
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

			//Parse variables in variable arrays
			for(var_xml = ezxml_child(var_arr_xml, "var"); var_xml; var_xml = var_xml->next){
				varname = ezxml_attr(var_xml, "name");
				varunits = ezxml_attr(var_xml, "units");
				vardesc = ezxml_attr(var_xml, "description");
				varstreams = ezxml_attr(var_xml, "streams");
				vararrgroup = ezxml_attr(var_xml, "array_group");
				varname_in_code = ezxml_attr(var_xml, "name_in_code");

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

				if(vararrpersistence == NULL){
					var_ptr->persistence = PERSISTENT;
				} else {
					if(strncmp(vararrpersistence, "persistent", 1024) == 0){
						var_ptr->persistence = PERSISTENT;
					} else if(strncmp(vararrpersistence, "scratch", 1024) == 0){
						var_ptr->persistence = SCRATCH;
					}
				}

				if(strncmp(vararrtype, "real", 1024) == 0){
					var_ptr->vtype = REAL;
				} else if(strncmp(vararrtype, "integer", 1024) == 0){
					var_ptr->vtype = INTEGER;
				} else if(strncmp(vararrtype, "logical", 1024) == 0){
					var_ptr->vtype = LOGICAL;
				} else if(strncmp(vararrtype, "text", 1024) == 0){
					var_ptr->vtype = CHARACTER;
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

				if(varstreams != NULL){
					snprintf(streams_buffer, 128, "%s", varstreams);
					if(strchr(streams_buffer, (int)'i')) var_ptr->iostreams |= INPUT0;
					if(strchr(streams_buffer, (int)'s')) var_ptr->iostreams |= SFC0;
					if(strchr(streams_buffer, (int)'r')) var_ptr->iostreams |= RESTART0;
					if(strchr(streams_buffer, (int)'o')) var_ptr->iostreams |= OUTPUT0;
				}

				if(varname_in_code == NULL){
					snprintf(var_ptr->name_in_code, 1024, "%s", varname);
				} else {
					snprintf(var_ptr->name_in_code, 1024, "%s", varname_in_code);
				}

				snprintf(var_ptr->super_array, 1024, "%s", vararrname);
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
			varstreams = ezxml_attr(var_xml, "streams");
			varname_in_code = ezxml_attr(var_xml, "name_in_code");

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

			if(varpersistence == NULL){
				var_ptr->persistence = PERSISTENT;
			} else {
				if(strncmp(varpersistence, "persistent", 1024) == 0){
					var_ptr->persistence = PERSISTENT;
				} else if(strncmp(varpersistence, "scratch", 1024) == 0){
					var_ptr->persistence = SCRATCH;
				}
			}

			if(strncmp(vartype, "real", 1024) == 0){
				var_ptr->vtype = REAL;
			} else if(strncmp(vartype, "integer", 1024) == 0){
				var_ptr->vtype = INTEGER;
			} else if(strncmp(vartype, "logical", 1024) == 0){
				var_ptr->vtype = LOGICAL;
			} else if(strncmp(vartype, "text", 1024) == 0){
				var_ptr->vtype = CHARACTER;
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

			if(varstreams != NULL){
				snprintf(streams_buffer, 128, "%s", varstreams);
				if(strchr(streams_buffer, (int)'i')) {
					var_ptr->iostreams |= INPUT0;
				}
				if(strchr(streams_buffer, (int)'s')) {
					var_ptr->iostreams |= SFC0;
				}
				if(strchr(streams_buffer, (int)'r')) {
					var_ptr->iostreams |= RESTART0;
				}
				if(strchr(streams_buffer, (int)'o')) {
					var_ptr->iostreams |= OUTPUT0;
				}
			}

			if(varname_in_code == NULL){
				snprintf(var_ptr->name_in_code, 1024, "%s", varname);
			} else {
				snprintf(var_ptr->name_in_code, 1024, "%s", varname_in_code);
			}

			snprintf(var_ptr->super_array, 1024, "-");
			snprintf(var_ptr->array_class, 1024, "-");

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

	return 0;
}

int getword(FILE * regfile, char * word)
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
}

int is_integer_constant(char * c) {
   int i;

   i = 0;
   while (c[i] != '\0') {
      if (c[i] < '0' || c[i] > '9') return -1;
      i++;
   }

   return atoi(c);
}

void sort_vars(struct variable * vars)
{
   struct variable * var_ptr;
   struct variable * var_ptr2;
   struct variable * var_ptr2_prev;
   char super_array[1024];
   char array_class[1024];

   var_ptr = vars;

/* Attempt at sorting first on super-array, then on class in the same loop
   while (var_ptr) {
      memcpy(super_array, var_ptr->super_array, 1024);
      memcpy(array_class, var_ptr->array_class, 1024);
      var_ptr2_prev = var_ptr;
      var_ptr2 = var_ptr->next;
      if (var_ptr2 && 
          (strncmp(super_array, var_ptr2->super_array, 1024) != 0 || strncmp(array_class, var_ptr2->array_class, 1024) != 0)) {
         while (var_ptr2) {
            if (strncmp(super_array, var_ptr2->super_array, 1024) == 0 && strncmp(array_class, var_ptr2->array_class, 1024) == 0) {
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
      memcpy(super_array, var_ptr->super_array, 1024);
      var_ptr2_prev = var_ptr;
      var_ptr2 = var_ptr->next;
      if (var_ptr2 && strncmp(super_array, var_ptr2->super_array, 1024) != 0) {
         while (var_ptr2) {
            if (strncmp(super_array, var_ptr2->super_array, 1024) == 0) {
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
}


void sort_group_vars(struct group_list * groups)
{
   struct variable_list * var_list;
   struct variable_list * var_ptr;
   struct variable_list * var_ptr2;
   struct variable_list * var_ptr2_prev;
   struct group_list * group_ptr;
   char super_array[1024];
   char array_class[1024];

   group_ptr = groups;

   while (group_ptr) {

      var_ptr = group_ptr->vlist;
   
      while (var_ptr) {
         memcpy(super_array, var_ptr->var->super_array, 1024);
         var_ptr2_prev = var_ptr;
         var_ptr2 = var_ptr->next;
         if (var_ptr2 != NULL && strncmp(super_array, var_ptr2->var->super_array, 1024) != 0) {
            while (var_ptr2) {
               if (strncmp(super_array, var_ptr2->var->super_array, 1024) == 0) {
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
}
