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
#include "dictionary.h"
#include "registry_types.h"
#include "gen_inc.h"
#include "fortprintf.h"

int is_derived_dim(char * d)
{
   if (strchr(d, (int)'+')) return 1;
   if (strchr(d, (int)'-')) return 1;

   return 0;
}


void get_outer_dim(struct variable * var, char * last_dim)
{
   struct dimension_list * dimlist_ptr;
 

   dimlist_ptr = var->dimlist;
   while (dimlist_ptr->next) dimlist_ptr = dimlist_ptr->next;

   strcpy(last_dim, dimlist_ptr->dim->name_in_file);
}

void split_derived_dim_string(char * dim, char ** p1, char ** p2)
{
   char * cp, * cm, * c;
   int n;

   cp = strchr(dim, (int)'+');
   cm = strchr(dim, (int)'-');
   if (!cp) 
      c = cm;
   else if (!cm) 
      c = cp;
   else if (cm < cp) 
      c = cm;
   else 
      c = cp;

   n = c - dim;
   *p1 = (char *)malloc(n*sizeof(char));
   snprintf(*p1, n, "%s", dim+1);

   *p2 = (char *)malloc((strlen(dim)-n+1)*sizeof(char));
   sprintf(*p2, "%s", dim+n);
}

void write_package_options(FILE *fd, struct package * pkgs){
	struct package * pkg_ptr;

	pkg_ptr = pkgs;

	fortprintf(fd, "%sActive", pkg_ptr->name);

	for (pkg_ptr = pkgs->next; pkg_ptr; pkg_ptr = pkg_ptr->next){
		fortprintf(fd, " .or. %sActive", pkg_ptr->name);
	}
}

void gen_namelists(struct namelist * nls)
{
   struct namelist * nls_ptr;
   struct dtable * dictionary;
   int done;
   char nlrecord[1024];
   FILE * fd;

   /*
    *  Generate config_type.inc
    */
   fd = fopen("config_defs.inc", "w");

   nls_ptr = nls;
   while (nls_ptr) {
      if (nls_ptr->vtype == INTEGER)   fortprintf(fd, "   integer :: %s\n",nls_ptr->name);
      if (nls_ptr->vtype == REAL)      fortprintf(fd, "   real (KIND=RKIND) :: %s\n",nls_ptr->name);
      if (nls_ptr->vtype == LOGICAL)   fortprintf(fd, "   logical :: %s\n",nls_ptr->name);
      if (nls_ptr->vtype == CHARACTER) fortprintf(fd, "   character (len=StrKIND) :: %s\n",nls_ptr->name);

      nls_ptr = nls_ptr->next;
   }

   fclose(fd);


   /*
    *  Generate namelist_defs.inc
    */
   fd = fopen("config_namelist_defs.inc", "w");
   dict_alloc(&dictionary);

   done = 0;
  
   while (!done) {
      nls_ptr = nls;
      while (nls_ptr && dict_search(dictionary, nls_ptr->record))
         nls_ptr = nls_ptr->next;

      if (nls_ptr) {
         dict_insert(dictionary, nls_ptr->record);
         strncpy(nlrecord, nls_ptr->record, 1024);
         fortprintf(fd, "      namelist /%s/ %s", nls_ptr->record, nls_ptr->name);
         nls_ptr = nls_ptr->next;
         while(nls_ptr) {
            if (strncmp(nls_ptr->record, nlrecord, 1024) == 0)
               fortprintf(fd, ", &\n                    %s", nls_ptr->name);
            nls_ptr = nls_ptr->next;
         }
         fortprintf(fd, "\n");
      }
      else
         done = 1;
   }
   

   dict_free(&dictionary);
   fclose(fd);


   /*
    *  Generate namelist_reads.inc
    */
   fd = fopen("config_set_defaults.inc", "w");
   nls_ptr = nls;
   while (nls_ptr) {
      if (nls_ptr->vtype == INTEGER) fortprintf(fd, "      %s = %i\n", nls_ptr->name, nls_ptr->defval.ival);
      if (nls_ptr->vtype == REAL)    fortprintf(fd, "      %s = %f\n", nls_ptr->name, nls_ptr->defval.rval);
      if (nls_ptr->vtype == LOGICAL) {
         if (nls_ptr->defval.lval == 0) 
            fortprintf(fd, "      %s = .false.\n", nls_ptr->name);
         else
            fortprintf(fd, "      %s = .true.\n", nls_ptr->name);
      }
      if (nls_ptr->vtype == CHARACTER)
         fortprintf(fd, "      %s = \"%s\"\n", nls_ptr->name, nls_ptr->defval.cval);
      nls_ptr = nls_ptr->next;
   }
   fortprintf(fd, "\n");
   fclose(fd);


   fd = fopen("config_namelist_reads.inc", "w");
   dict_alloc(&dictionary);
   nls_ptr = nls;
   while (nls_ptr) {
      if (!dict_search(dictionary, nls_ptr->record)) {
         fortprintf(fd, "         read(funit,%s,iostat=ierr)\n", nls_ptr->record);
         fortprintf(fd, "         if (ierr > 0) then\n");
         fortprintf(fd, "            write(0,*) \'Error while reading namelist record &%s\'\n",nls_ptr->record);
         fortprintf(fd, "            call mpas_dmpar_abort(dminfo)\n");
         fortprintf(fd, "         else if (ierr < 0) then\n");
         fortprintf(fd, "            write(0,*) \'Namelist record &%s not found; using default values for this namelist\'\'s variables\'\n",nls_ptr->record);
         fortprintf(fd, "         end if\n");
         fortprintf(fd, "         rewind(funit)\n");

         dict_insert(dictionary, nls_ptr->record);
      }
      nls_ptr = nls_ptr->next;
   }
   fortprintf(fd, "\n");

   dict_free(&dictionary);
   fclose(fd);


   fd = fopen("config_bcast_namelist.inc", "w");
   nls_ptr = nls;
   while (nls_ptr) {
      if (nls_ptr->vtype == INTEGER)   fortprintf(fd, "      call mpas_dmpar_bcast_int(dminfo, %s)\n", nls_ptr->name);
      if (nls_ptr->vtype == REAL)      fortprintf(fd, "      call mpas_dmpar_bcast_real(dminfo, %s)\n", nls_ptr->name);
      if (nls_ptr->vtype == LOGICAL)   fortprintf(fd, "      call mpas_dmpar_bcast_logical(dminfo, %s)\n", nls_ptr->name);
      if (nls_ptr->vtype == CHARACTER) fortprintf(fd, "      call mpas_dmpar_bcast_char(dminfo, %s)\n", nls_ptr->name);
      nls_ptr = nls_ptr->next;
   }
   fortprintf(fd, "\n");
   fclose(fd);
}

void gen_history_attributes(char * modelname, char * corename, char * version)
{
	FILE *fd;

	fd = fopen("model_variables.inc","w");
	fortprintf(fd, "       character (len=StrKIND) :: modelName = '%s' !< Constant: Name of model\n", modelname);
	fortprintf(fd, "       character (len=StrKIND) :: coreName = '%s' !< Constant: Name of core\n", corename);
	fortprintf(fd, "       character (len=StrKIND) :: modelVersion = '%s' !< Constant: Version number\n", version);
	fclose(fd);
}


void gen_field_defs(struct group_list * groups, struct variable * vars, struct dimension * dims)
{
   struct variable * var_ptr;
   struct variable * var_ptr2;
   struct variable_list * var_list_ptr;
   struct variable_list * var_list_ptr2;
   struct variable_list * var_list_ptr3;
   struct dimension * dim_ptr;
   struct dimension_list * dimlist_ptr;
   struct group_list * group_ptr;
   FILE * fd, *fd2;
   char var_array[1024];
   char array_class[1024];
   char outer_dim[1024];
   int i, new_class;
   int class_start, class_end;
   int vtype;
   char type_str[7];


   /*
    *  Generate declarations of dimensions
    */
   fd = fopen("field_dimensions.inc", "w");
   dim_ptr = dims;
   while (dim_ptr) {
      if (dim_ptr->constant_value < 0 && !dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) fortprintf(fd, "      integer :: %s\n", dim_ptr->name_in_code);
      if (dim_ptr->constant_value < 0 && dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) fortprintf(fd, "      integer :: %s\n", dim_ptr->name_in_file);
      dim_ptr = dim_ptr->next;
   }
   dim_ptr = dims;
   while (dim_ptr) {
      if (dim_ptr->constant_value < 0 && !dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) {
		  fortprintf(fd, "      integer :: %sSolve\n", dim_ptr->name_in_code);
		  fortprintf(fd, "      integer, dimension(:), pointer :: %sArray\n", dim_ptr->name_in_code);
	  }
      if (dim_ptr->constant_value < 0 && dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) {
		  fortprintf(fd, "      integer :: %sSolve\n", dim_ptr->name_in_file);
	  }
      dim_ptr = dim_ptr->next;
   }

   fclose(fd);

   /*
    *  Generate dummy dimension argument list
    */
   fd = fopen("dim_dummy_args.inc", "w");
   dim_ptr = dims;
   if (dim_ptr && dim_ptr->constant_value < 0 && !dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) {
      fortprintf(fd, "                            %s", dim_ptr->name_in_code);
      dim_ptr = dim_ptr->next;
   }
   else if (dim_ptr && dim_ptr->constant_value < 0 && dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) {
      fortprintf(fd, "                            %s", dim_ptr->name_in_file);
      dim_ptr = dim_ptr->next;
   }
   while (dim_ptr) {
      if (dim_ptr->constant_value < 0 && !dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) fortprintf(fd, ", %s", dim_ptr->name_in_code);
      if (dim_ptr->constant_value < 0 && dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) fortprintf(fd, ", %s", dim_ptr->name_in_file);
      dim_ptr = dim_ptr->next;
   }
   fortprintf(fd, " &\n");

   fclose(fd);

   /*
    *  Generate dummy dimension argument declaration list
    */
   fd = fopen("dim_dummy_decls.inc", "w");
   dim_ptr = dims;
   if (dim_ptr && dim_ptr->constant_value < 0 && !dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) {
      fortprintf(fd, "      integer, intent(in) :: %s", dim_ptr->name_in_code);
      dim_ptr = dim_ptr->next;
   }
   else if (dim_ptr && dim_ptr->constant_value < 0 && dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) {
      fortprintf(fd, "      integer, intent(in) :: %s", dim_ptr->name_in_file);
      dim_ptr = dim_ptr->next;
   }
   while (dim_ptr) {
      if (dim_ptr->constant_value < 0 && !dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) fortprintf(fd, ", %s", dim_ptr->name_in_code);
      if (dim_ptr->constant_value < 0 && dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) fortprintf(fd, ", %s", dim_ptr->name_in_file);
      dim_ptr = dim_ptr->next;
   }
   fortprintf(fd, "\n");

   fclose(fd);

   /*
    *  Generate dummy dimension argument declaration list
    */
   fd = fopen("dim_dummy_decls_inout.inc", "w");
   dim_ptr = dims;
   if (dim_ptr && dim_ptr->constant_value < 0 && !dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) {
      fortprintf(fd, "      integer, intent(inout) :: %s", dim_ptr->name_in_code);
      dim_ptr = dim_ptr->next;
   }
   else if (dim_ptr && dim_ptr->constant_value < 0 && dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) {
      fortprintf(fd, "      integer, intent(inout) :: %s", dim_ptr->name_in_file);
      dim_ptr = dim_ptr->next;
   }
   while (dim_ptr) {
      if (dim_ptr->constant_value < 0 && !dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) fortprintf(fd, ", %s", dim_ptr->name_in_code);
      if (dim_ptr->constant_value < 0 && dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) fortprintf(fd, ", %s", dim_ptr->name_in_file);
      dim_ptr = dim_ptr->next;
   }
   fortprintf(fd, "\n");

   fclose(fd);

   /*
    *  Generate non-input dummy dimension argument declaration list
    */
   fd = fopen("dim_dummy_decls_noinput.inc", "w");
   dim_ptr = dims;
   if (dim_ptr && dim_ptr->constant_value < 0 && !dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) {
      fortprintf(fd, "      integer :: %s", dim_ptr->name_in_code);
      dim_ptr = dim_ptr->next;
   }
   else if (dim_ptr && dim_ptr->constant_value < 0 && dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) {
      fortprintf(fd, "      integer :: %s", dim_ptr->name_in_file);
      dim_ptr = dim_ptr->next;
   }
   while (dim_ptr) {
      if (dim_ptr->constant_value < 0 && !dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) fortprintf(fd, ", %s", dim_ptr->name_in_code);
      if (dim_ptr->constant_value < 0 && dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) fortprintf(fd, ", %s", dim_ptr->name_in_file);
      dim_ptr = dim_ptr->next;
   }
   fortprintf(fd, "\n");

   fclose(fd);



   /*
    *  Generate dummy dimension assignment instructions
    */
   fd = fopen("dim_dummy_assigns.inc", "w");
   dim_ptr = dims;
   if (dim_ptr && dim_ptr->constant_value < 0 && !dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) {
      fortprintf(fd, "      %s = block %% mesh %% %s\n", dim_ptr->name_in_code, dim_ptr->name_in_code);
      dim_ptr = dim_ptr->next;
   } 
   else if (dim_ptr && dim_ptr->constant_value < 0 && dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) {
      fortprintf(fd, "      %s = block %% mesh %% %s\n", dim_ptr->name_in_file, dim_ptr->name_in_file);
      dim_ptr = dim_ptr->next;
   }
   while (dim_ptr) {
      if (dim_ptr->constant_value < 0 && !dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) fortprintf(fd, "      %s = block %% mesh %% %s\n", dim_ptr->name_in_code, dim_ptr->name_in_code);
      if (dim_ptr->constant_value < 0 && dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) fortprintf(fd, "      %s = block %% mesh %% %s\n", dim_ptr->name_in_file, dim_ptr->name_in_file);
      dim_ptr = dim_ptr->next;
   }
   fortprintf(fd, "\n");

   fclose(fd);


   /*
    *  Generate declarations of dimensions
    */
   fd = fopen("dim_decls.inc", "w");
   dim_ptr = dims;
   while (dim_ptr) {
      if (dim_ptr->constant_value < 0 && !dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) fortprintf(fd, "      integer :: %s\n", dim_ptr->name_in_code);
      if (dim_ptr->constant_value < 0 && dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) fortprintf(fd, "      integer :: %s\n", dim_ptr->name_in_file);
      dim_ptr = dim_ptr->next;
   }

   fclose(fd);


   /*
    *  Generate calls to read dimensions from input file
    */
   fd = fopen("read_dims.inc", "w");
   dim_ptr = dims;
   while (dim_ptr) {
      if (dim_ptr->constant_value < 0 && !dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) fortprintf(fd, "      call MPAS_io_inq_dim(inputHandle, \'%s\', %s, ierr)\n", dim_ptr->name_in_file, dim_ptr->name_in_code);
      else if (dim_ptr->constant_value < 0 && dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) fortprintf(fd, "      %s = %s\n", dim_ptr->name_in_file, dim_ptr->name_in_code);
      dim_ptr = dim_ptr->next;
   }

   fclose(fd);


   /*
    *  Generate declarations of mesh group
    */
   fd = fopen("time_invariant_fields.inc", "w");
   group_ptr = groups;
   while (group_ptr) {

      if (!strncmp(group_ptr->name, "mesh", 1024)) {

         var_list_ptr = group_ptr->vlist;
         memcpy(var_array, var_list_ptr->var->var_array, 1024);
         i = 1;
         while (var_list_ptr) {
            if (strncmp(var_array, var_list_ptr->var->var_array, 1024) != 0) {
               memcpy(var_array, var_list_ptr->var->var_array, 1024);
               i = 1;
             }
            if (strncmp(var_list_ptr->var->array_class, "-", 1024) != 0) fortprintf(fd, "      integer :: index_%s = -1\n", var_list_ptr->var->name_in_code);
            var_list_ptr = var_list_ptr->next;
         }

         var_list_ptr = group_ptr->vlist;
         memcpy(var_array, var_list_ptr->var->var_array, 1024);
         memcpy(array_class, var_list_ptr->var->array_class, 1024);
         class_start = 1;
         class_end = 1;
         i = 1;
         while (var_list_ptr) {
            if (strncmp(var_list_ptr->var->var_array, "-", 1024) != 0) {
               if (strncmp(var_array, var_list_ptr->var->var_array, 1024) != 0) {
                  if (strncmp(var_array, "-", 1024) != 0) fortprintf(fd, "      integer :: %s_end = -1\n", array_class);
                  if (strncmp(var_array, "-", 1024) != 0) fortprintf(fd, "      integer :: num_%s = -1\n", var_array);
                  class_start = 1;
                  class_end = 1;
                  i = 1;
                  memcpy(var_array, var_list_ptr->var->var_array, 1024);
                  memcpy(array_class, var_list_ptr->var->array_class, 1024);
                  fortprintf(fd, "      integer :: %s_start = -1\n", array_class);
               }
               else if (strncmp(array_class, var_list_ptr->var->array_class, 1024) != 0) {
                  fortprintf(fd, "      integer :: %s_end = -1\n", array_class);
                  class_start = class_end+1;
                  class_end = class_start;
                  memcpy(array_class, var_list_ptr->var->array_class, 1024);
                  fortprintf(fd, "      integer :: %s_start = -1\n", array_class);
                  i++;
               }
               else {
                  class_end++;
                  i++;
               }
            }
            var_list_ptr = var_list_ptr->next;
         }
         if (strncmp(var_array, "-", 1024) != 0) fortprintf(fd, "      integer :: %s_end = -1\n", array_class);
         if (strncmp(var_array, "-", 1024) != 0) fortprintf(fd, "      integer :: num_%s = -1\n", var_array);

         var_list_ptr = group_ptr->vlist;
         while (var_list_ptr) {
            var_ptr = var_list_ptr->var;
            if (!strncmp(var_ptr->var_array, "-", 1024)) {
              if (var_ptr->vtype == INTEGER)   fortprintf(fd, "      type (field%idInteger), pointer :: %s\n", var_ptr->ndims, var_ptr->name_in_code);
              if (var_ptr->vtype == REAL)      fortprintf(fd, "      type (field%idReal), pointer :: %s\n", var_ptr->ndims, var_ptr->name_in_code);
              if (var_ptr->vtype == CHARACTER) fortprintf(fd, "      type (field%idChar), pointer :: %s\n", var_ptr->ndims, var_ptr->name_in_code);
            }
            else {
              if (var_ptr->vtype == INTEGER)   fortprintf(fd, "      type (field%idInteger), pointer :: %s\n", var_ptr->ndims+1, var_ptr->var_array);
              if (var_ptr->vtype == REAL)      fortprintf(fd, "      type (field%idReal), pointer :: %s\n", var_ptr->ndims+1, var_ptr->var_array);
              if (var_ptr->vtype == CHARACTER) fortprintf(fd, "      type (field%idChar), pointer :: %s\n", var_ptr->ndims+1, var_ptr->var_array);
              while (var_list_ptr->next && !strncmp(var_list_ptr->next->var->var_array, var_list_ptr->var->var_array, 1024)) var_list_ptr = var_list_ptr->next;
            }
            var_list_ptr = var_list_ptr->next;
         }
         break;
      }
      group_ptr = group_ptr->next;
   }

   fclose(fd);


   /*
    *  Generate declarations of non-mesh groups
    */
   fd = fopen("variable_groups.inc", "w");

   group_ptr = groups;
   while (group_ptr) {
      if (strncmp(group_ptr->name, "mesh", 1024)) {
         fortprintf(fd, "   type %s_type\n", group_ptr->name);

         fortprintf(fd, "      type (block_type), pointer :: block\n");

         var_list_ptr = group_ptr->vlist;
         if(group_ptr->vlist != NULL){
            memcpy(var_array, var_list_ptr->var->var_array, 1024);
            i = 1;
            while (var_list_ptr) {
               if (strncmp(var_array, var_list_ptr->var->var_array, 1024) != 0) {
                  memcpy(var_array, var_list_ptr->var->var_array, 1024);
                  i = 1;
                }
               if (strncmp(var_list_ptr->var->array_class, "-", 1024) != 0) fortprintf(fd, "      integer :: index_%s = -1\n", var_list_ptr->var->name_in_code);
               var_list_ptr = var_list_ptr->next;
            }

            var_list_ptr = group_ptr->vlist;
            sprintf(var_array, "-");
            sprintf(array_class, "-");
            class_start = 1;
            class_end = 1;
            i = 1;

            while (var_list_ptr) {

               /* Is the current variable in a super array? */
               if (strncmp(var_list_ptr->var->var_array, "-", 1024) != 0) {

                  /* Have we hit the beginning of a new super array? */
                  if (strncmp(var_array, var_list_ptr->var->var_array, 1024) != 0) {
                     /* Finish off the previous super array? */
                     if (strncmp(var_array, "-", 1024) != 0) {
                        fortprintf(fd, "      integer :: %s_end = -1\n", array_class);
                        fortprintf(fd, "      integer :: num_%s = -1\n", var_array);
                     }
                     class_start = 1;
                     class_end = 1;
                     i = 1;
                     memcpy(var_array, var_list_ptr->var->var_array, 1024);
                     memcpy(array_class, var_list_ptr->var->array_class, 1024);
                     fortprintf(fd, "      integer :: %s_start = -1\n", array_class);
                  }
                  /* Or have we hit the beginning of a new array class? */
                  else if (strncmp(array_class, var_list_ptr->var->array_class, 1024) != 0) {
                     fortprintf(fd, "      integer :: %s_end = -1\n", array_class);
                     class_start = class_end+1;
                     class_end = class_start;
                     memcpy(array_class, var_list_ptr->var->array_class, 1024);
                     fortprintf(fd, "      integer :: %s_start = -1\n", array_class);
                     i++;
                  }
                  else {
                     class_end++;
                     i++;
                  }

               }
               var_list_ptr = var_list_ptr->next;

            }
            if (strncmp(var_array, "-", 1024) != 0) fortprintf(fd, "      integer :: %s_end = -1\n", array_class);
            if (strncmp(var_array, "-", 1024) != 0) fortprintf(fd, "      integer :: num_%s = -1\n", var_array);

            var_list_ptr = group_ptr->vlist;
            while (var_list_ptr) {
               var_ptr = var_list_ptr->var;
               if (!strncmp(var_ptr->var_array, "-", 1024)) {
                 if (var_ptr->vtype == INTEGER)   fortprintf(fd, "      type (field%idInteger), pointer :: %s\n", var_ptr->ndims, var_ptr->name_in_code);
                 if (var_ptr->vtype == REAL)      fortprintf(fd, "      type (field%idReal), pointer :: %s\n", var_ptr->ndims, var_ptr->name_in_code);
                 if (var_ptr->vtype == CHARACTER) fortprintf(fd, "      type (field%idChar), pointer :: %s\n", var_ptr->ndims, var_ptr->name_in_code);
               }
               else {
                 if (var_ptr->vtype == INTEGER)   fortprintf(fd, "      type (field%idInteger), pointer :: %s\n", var_ptr->ndims+1, var_ptr->var_array);
                 if (var_ptr->vtype == REAL)      fortprintf(fd, "      type (field%idReal), pointer :: %s\n", var_ptr->ndims+1, var_ptr->var_array);
                 if (var_ptr->vtype == CHARACTER) fortprintf(fd, "      type (field%idChar), pointer :: %s\n", var_ptr->ndims+1, var_ptr->var_array);
                 while (var_list_ptr->next && !strncmp(var_list_ptr->next->var->var_array, var_list_ptr->var->var_array, 1024)) var_list_ptr = var_list_ptr->next;
               }
               var_list_ptr = var_list_ptr->next;
            }
         }
   
         fortprintf(fd, "   end type %s_type\n\n\n", group_ptr->name);
   
         if (group_ptr->ntime_levs > 1) {
            fortprintf(fd, "   type %s_pointer_type\n", group_ptr->name);
            fortprintf(fd, "      type (%s_type), pointer :: %s \n", group_ptr->name, group_ptr->name);
            fortprintf(fd, "   end type %s_pointer_type\n\n\n", group_ptr->name);
   
            fortprintf(fd, "   type %s_multilevel_type\n", group_ptr->name);
            fortprintf(fd, "      integer :: nTimeLevels\n");
            fortprintf(fd, "      type (%s_pointer_type), dimension(:), pointer :: time_levs\n", group_ptr->name);
            fortprintf(fd, "   end type %s_multilevel_type\n\n\n", group_ptr->name);
         }
      }
      group_ptr = group_ptr->next;
   }

   fclose(fd);

   /*
    *  Generate instantiations of variable groups in block_type
    */
   fd = fopen("block_group_members.inc", "w");

   group_ptr = groups;
   while (group_ptr) {
      if (group_ptr->ntime_levs > 1) {
         fortprintf(fd, "      type (%s_multilevel_type), pointer :: %s\n", group_ptr->name, group_ptr->name);
         fortprintf(fd, "      type (%s_type), pointer :: provis_%s\n", group_ptr->name, group_ptr->name);
	  } else {
         fortprintf(fd, "      type (%s_type), pointer :: %s\n", group_ptr->name, group_ptr->name);
	  }
      group_ptr = group_ptr->next;
   }

   fclose(fd);


   /*
    *  Generate routines for allocating provisional types
    */
   fd = fopen("provis_alloc_routines.inc", "w");

   group_ptr = groups;
   while (group_ptr) {
      if (group_ptr->ntime_levs > 1) {
	 fortprintf(fd, "   subroutine mpas_setup_provis_%s(b)!{{{\n", group_ptr->name);
	 fortprintf(fd, "      type (block_type), pointer :: b\n");
	 fortprintf(fd, "      type (block_type), pointer :: block\n\n");
	 fortprintf(fd, "#include \"dim_dummy_decls_noinput.inc\"\n\n");
	 fortprintf(fd, "      block => b\n");
	 fortprintf(fd, "      do while(associated(block))\n");
	 fortprintf(fd, "#include \"dim_dummy_assigns.inc\"\n\n");
         fortprintf(fd, "         allocate(block %% provis_%s)\n", group_ptr->name);
         fortprintf(fd, "         call mpas_allocate_%s(block, block %% provis_%s, &\n", group_ptr->name, group_ptr->name);
         fortprintf(fd, "#include \"dim_dummy_args.inc\"\n");
         fortprintf(fd, "                              )\n\n");
	 fortprintf(fd, "         block => block %% next \n");
	 fortprintf(fd, "      end do\n\n");
	 fortprintf(fd, "      block => b\n");
	 fortprintf(fd, "      do while(associated(block))\n");
         fortprintf(fd, "         if(associated(block %% prev) .and. associated(block %% next)) then\n");
         fortprintf(fd, "            call mpas_create_%s_links(block %% provis_%s, prev = block %% prev %% provis_%s, next = block %% next %% provis_%s)\n", group_ptr->name, group_ptr->name, group_ptr->name, group_ptr->name);
         fortprintf(fd, "         else if(associated(block %% prev)) then\n");
         fortprintf(fd, "            call mpas_create_%s_links(block %% provis_%s, prev = block %% prev %% provis_%s)\n", group_ptr->name, group_ptr->name, group_ptr->name);
         fortprintf(fd, "         else if(associated(block %% next)) then\n");
         fortprintf(fd, "            call mpas_create_%s_links(block %% provis_%s, next = block %% next %% provis_%s)\n", group_ptr->name, group_ptr->name, group_ptr->name);
         fortprintf(fd, "         else\n");
         fortprintf(fd, "            call mpas_create_%s_links(block %% provis_%s)\n", group_ptr->name, group_ptr->name);
         fortprintf(fd, "         end if\n");
	 fortprintf(fd, "         block => block %% next \n");
	 fortprintf(fd, "      end do\n");
	 fortprintf(fd, "   end subroutine mpas_setup_provis_%s!}}}\n\n", group_ptr->name);

	 fortprintf(fd, "   subroutine mpas_deallocate_provis_%s(b)!{{{\n", group_ptr->name);
	 fortprintf(fd, "      type (block_type), pointer :: b\n");
	 fortprintf(fd, "      type (block_type), pointer :: block\n\n");
	 fortprintf(fd, "      block => b\n");
	 fortprintf(fd, "      do while(associated(block))\n");
	 fortprintf(fd, "         call mpas_deallocate_%s(block %% provis_%s)\n", group_ptr->name, group_ptr->name);
	 fortprintf(fd, "         deallocate(block %% provis_%s)\n", group_ptr->name);
	 fortprintf(fd, "         block => block %% next\n");
	 fortprintf(fd, "      end do\n");
	 fortprintf(fd, "   end subroutine mpas_deallocate_provis_%s!}}}\n", group_ptr->name);
	  }
      group_ptr = group_ptr->next;
   }
   fclose(fd);



   /* To be included in allocate_block */
   fd = fopen("block_allocs.inc", "w");
   group_ptr = groups;
   while (group_ptr) {
      fortprintf(fd, "      allocate(b %% %s)\n", group_ptr->name);
      if (group_ptr->ntime_levs > 1) {
         fortprintf(fd, "      b %% %s %% nTimeLevels = %i\n", group_ptr->name, group_ptr->ntime_levs);
         fortprintf(fd, "      allocate(b %% %s %% time_levs(%i))\n", group_ptr->name, group_ptr->ntime_levs);
         fortprintf(fd, "      do i=1,b %% %s %% nTimeLevels\n", group_ptr->name);
         fortprintf(fd, "         allocate(b %% %s %% time_levs(i) %% %s)\n", group_ptr->name, group_ptr->name);
         fortprintf(fd, "         call mpas_allocate_%s(b, b %% %s %% time_levs(i) %% %s, &\n", group_ptr->name, group_ptr->name, group_ptr->name);
         fortprintf(fd, "#include \"dim_dummy_args.inc\"\n");
         fortprintf(fd, "                         )\n");
         fortprintf(fd, "      end do\n\n");
      }
      else {
         fortprintf(fd, "      call mpas_allocate_%s(b, b %% %s, &\n", group_ptr->name, group_ptr->name);
         fortprintf(fd, "#include \"dim_dummy_args.inc\"\n");
         fortprintf(fd, "                      )\n\n");
      }
      group_ptr = group_ptr->next;
   }
   fclose(fd);

   
   /* To be included in deallocate_block */
   fd = fopen("block_deallocs.inc", "w");
   group_ptr = groups;
   while (group_ptr) {
      if (group_ptr->ntime_levs > 1) {
         fortprintf(fd, "      do i=1,b %% %s %% nTimeLevels\n", group_ptr->name);
         fortprintf(fd, "         call mpas_deallocate_%s(b %% %s %% time_levs(i) %% %s)\n", group_ptr->name, group_ptr->name, group_ptr->name);
         fortprintf(fd, "         deallocate(b %% %s %% time_levs(i) %% %s)\n", group_ptr->name, group_ptr->name);
         fortprintf(fd, "      end do\n");
         fortprintf(fd, "      deallocate(b %% %s %% time_levs)\n", group_ptr->name);
      }
      else {
         fortprintf(fd, "      call mpas_deallocate_%s(b %% %s)\n", group_ptr->name, group_ptr->name);
      }
      fortprintf(fd, "      deallocate(b %% %s)\n\n", group_ptr->name);
      group_ptr = group_ptr->next;
   }
   fclose(fd);

   /* Definitions of allocate subroutines */
   fd = fopen("group_alloc_routines.inc", "w");
   group_ptr = groups;
   while (group_ptr) {
      fortprintf(fd, "   subroutine mpas_allocate_%s(b, %s, &\n", group_ptr->name, group_ptr->name);
      fortprintf(fd, "#include \"dim_dummy_args.inc\"\n");
      fortprintf(fd, "                         )\n");
      fortprintf(fd, "\n");
      fortprintf(fd, "      implicit none\n");
      fortprintf(fd, "\n");
      fortprintf(fd, "      type (block_type), pointer :: b\n");
      fortprintf(fd, "      type (%s_type), intent(inout) :: %s\n", group_ptr->name, group_ptr->name);
      fortprintf(fd, "      logical :: group_started\n");
      fortprintf(fd, "      integer :: index_counter\n");
      fortprintf(fd, "      integer :: group_counter\n");
      fortprintf(fd, "      integer :: group_start\n");
      fortprintf(fd, "#include \"dim_dummy_decls.inc\"\n");
      fortprintf(fd, "\n");

      fortprintf(fd, "      %s %% block => b\n", group_ptr->name);

      if (!strncmp(group_ptr->name, "mesh", 1024)) {
         dim_ptr = dims;
         while (dim_ptr) {
            if (dim_ptr->constant_value < 0 && !dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) fortprintf(fd, "      %s %% %s = %s\n", group_ptr->name, dim_ptr->name_in_code, dim_ptr->name_in_code);
            if (dim_ptr->constant_value < 0 && dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) fortprintf(fd, "      %s %% %s = %s\n", group_ptr->name, dim_ptr->name_in_file, dim_ptr->name_in_file);
            dim_ptr = dim_ptr->next;
         }

         fortprintf(fd, "\n");
      }


      var_list_ptr = group_ptr->vlist;
      while (var_list_ptr) {
         var_ptr = var_list_ptr->var;
         if (strncmp(var_ptr->var_array, "-", 1024) != 0) {
            memcpy(var_array, var_ptr->var_array, 1024);
            memcpy(array_class, var_ptr->array_class, 1024);
            i = 0;
            while (var_list_ptr && strncmp(var_array, var_list_ptr->var->var_array, 1024) == 0) {
               i++;
               var_list_ptr2 = var_list_ptr;
               var_list_ptr = var_list_ptr->next;
            }
            var_ptr2 = var_list_ptr2->var;
            fortprintf(fd, "      index_counter = 0\n");
            fortprintf(fd, "      group_counter = -1\n");
            fortprintf(fd, "      group_start = -1\n");
            fortprintf(fd, "      group_started = .false.\n");
            fortprintf(fd, "      allocate(%s %% %s)\n", group_ptr->name, var_ptr2->var_array);
            fortprintf(fd, "      allocate(%s %% %s %% ioinfo)\n", group_ptr->name, var_ptr2->var_array);
            fortprintf(fd, "      %s %% %s %% fieldName = \'%s\'\n", group_ptr->name, var_ptr2->var_array, var_ptr2->var_array);
            fortprintf(fd, "      %s %% %s %% isVarArray = .true.\n", group_ptr->name, var_ptr2->var_array);
			/* Initialization of indices and size */
			i = 0;
			new_class = 0;
            var_list_ptr3 = group_ptr->vlist;
            memcpy(array_class, "-", 1024);
            while (var_list_ptr3) {
               if (strncmp(var_array, var_list_ptr3->var->var_array, 1024) == 0) {
			      if (strncmp(array_class, var_list_ptr3->var->array_class, 1024) != 0) {
					  if (strncmp(array_class, "-", 1024) != 0) {
						  fortprintf(fd, "      if (group_counter >= 0) then\n");
						  fortprintf(fd, "         %s %% %s_start = group_start\n", group_ptr->name, array_class);
						  fortprintf(fd, "         %s %% %s_end = group_start + group_counter\n", group_ptr->name, array_class);
						  fortprintf(fd, "      end if\n");
						  fortprintf(fd, "      group_counter = -1\n");
						  fortprintf(fd, "      group_started = .false.\n");
					  }
					  memcpy(array_class, var_list_ptr3->var->array_class, 1024);
				  }

			      if(var_list_ptr3->var->persistence != PACKAGE){
					  fortprintf(fd, "      index_counter = index_counter + 1\n");
					  fortprintf(fd, "      group_counter = group_counter + 1\n");
					  fortprintf(fd, "      %s %% index_%s = index_counter\n", group_ptr->name, var_list_ptr3->var->name_in_code);
					  fortprintf(fd, "      if (.not. group_started) then\n");
					  fortprintf(fd, "         group_start = index_counter\n");
					  fortprintf(fd, "         group_started = .true.\n");
					  fortprintf(fd, "      end if\n");
				  } else {
					  fortprintf(fd, "      if (");
					  write_package_options(fd, var_list_ptr3->var->package_name);
					  fortprintf(fd, ") then\n");

					  fortprintf(fd, "         index_counter = index_counter + 1\n");
					  fortprintf(fd, "         group_counter = group_counter + 1\n");
					  fortprintf(fd, "         %s %% index_%s = index_counter\n", group_ptr->name, var_list_ptr3->var->name_in_code);
					  fortprintf(fd, "         if (.not. group_started) then\n");
					  fortprintf(fd, "            group_start = index_counter\n");
					  fortprintf(fd, "            group_started = .true.\n");
					  fortprintf(fd, "         end if\n");
					  fortprintf(fd, "      end if\n");
				  }
               }
               var_list_ptr3 = var_list_ptr3->next;
            }

			fortprintf(fd, "      if (group_counter > 0) then\n");
			fortprintf(fd, "         %s %% %s_start = group_start\n", group_ptr->name, array_class);
			fortprintf(fd, "         %s %% %s_end = group_start + group_counter\n", group_ptr->name, array_class);
			fortprintf(fd, "      end if\n");
			fortprintf(fd, "      %s %% num_%s = index_counter\n", group_ptr->name, var_array);
			fortprintf(fd, "      if ( %s %% num_%s > 0 ) then\n", group_ptr->name, var_array);
			fortprintf(fd, "         allocate(%s %% %s %% constituentNames(%s %% num_%s))\n", group_ptr->name, var_array, group_ptr->name, var_array);
			fortprintf(fd, "      end if\n");

            /* Initialization for constituent names */
            i = 0;
            var_list_ptr3 = group_ptr->vlist;
            while (var_list_ptr3) {
               if (strncmp(var_array, var_list_ptr3->var->var_array, 1024) == 0) {
			      if(var_list_ptr3->var->persistence != PACKAGE) {
					  fortprintf(fd, "      %s %% %s %% constituentNames(%s %% index_%s) = \'%s\'\n", group_ptr->name, var_array, group_ptr->name, var_list_ptr3->var->name_in_code, var_list_ptr3->var->name_in_file);
				  } else {
					  fortprintf(fd, "      if (%s %% index_%s > 0) then\n", group_ptr->name, var_list_ptr3->var->name_in_code);
					  fortprintf(fd, "         %s %% %s %% constituentNames(%s %% index_%s) = \'%s\'\n", group_ptr->name, var_array, group_ptr->name, var_list_ptr3->var->name_in_code, var_list_ptr3->var->name_in_file);
					  fortprintf(fd, "      end if\n");

				  }
               }
               var_list_ptr3 = var_list_ptr3->next;
            }

			if(var_ptr2->persistence == PERSISTENT || var_ptr2->persistence == PACKAGE){
			   fortprintf(fd, "      %s %% %s %% isPersistent = .true.\n", group_ptr->name, var_ptr2->var_array);
			   if(var_ptr2->persistence == PACKAGE){
                 fortprintf(fd, "      if (%s %% num_%s > 0) then\n", group_ptr->name, var_ptr2->var_array);
				 fortprintf(fd, "         %s %% %s %% isActive = .true.\n", group_ptr->name, var_ptr2->var_array);
				 fortprintf(fd, "         allocate(%s %% %s %% array(%s %% num_%s, ", group_ptr->name, var_ptr2->var_array, group_ptr->name, var_ptr2->var_array);
			   } else {
				 fortprintf(fd, "      %s %% %s %% isActive = .true.\n", group_ptr->name, var_ptr2->var_array);
				 fortprintf(fd, "      allocate(%s %% %s %% array(%s %% num_%s, ", group_ptr->name, var_ptr2->var_array, group_ptr->name, var_ptr2->var_array);
			   }
               dimlist_ptr = var_ptr2->dimlist;
               if (!strncmp(dimlist_ptr->dim->name_in_file, "nCells", 1024) ||
                   !strncmp(dimlist_ptr->dim->name_in_file, "nEdges", 1024) ||
                   !strncmp(dimlist_ptr->dim->name_in_file, "nVertices", 1024))
                  if (!dimlist_ptr->dim->namelist_defined) fortprintf(fd, "%s + 1", dimlist_ptr->dim->name_in_code);
                  else fortprintf(fd, "%s + 1", dimlist_ptr->dim->name_in_file);
               else
                  if (dimlist_ptr->dim->namelist_defined) fortprintf(fd, "%s", dimlist_ptr->dim->name_in_file);
                  else fortprintf(fd, "%s", dimlist_ptr->dim->name_in_code);
               dimlist_ptr = dimlist_ptr->next;
               while (dimlist_ptr) {
                  if (!strncmp(dimlist_ptr->dim->name_in_file, "nCells", 1024) ||
                      !strncmp(dimlist_ptr->dim->name_in_file, "nEdges", 1024) ||
                      !strncmp(dimlist_ptr->dim->name_in_file, "nVertices", 1024))
                     if (!dimlist_ptr->dim->namelist_defined) fortprintf(fd, ", %s + 1", dimlist_ptr->dim->name_in_code);
                     else fortprintf(fd, ", %s + 1", dimlist_ptr->dim->name_in_file);
                  else
                     if (dimlist_ptr->dim->namelist_defined) fortprintf(fd, ", %s", dimlist_ptr->dim->name_in_file);
                     else fortprintf(fd, ", %s", dimlist_ptr->dim->name_in_code);
                  dimlist_ptr = dimlist_ptr->next;
               }
               fortprintf(fd, "))\n");
			   if (var_ptr2->persistence == PACKAGE) {
				  fortprintf(fd, "         %s %% %s %% array = %s\n", group_ptr->name, var_ptr2->var_array, var_ptr2->default_value ); /* initialize field */
                  fortprintf(fd, "      else\n");
			      fortprintf(fd, "         %s %% %s %% isActive = .false.\n", group_ptr->name, var_ptr2->var_array);
                  fortprintf(fd, "         nullify(%s %% %s %% array)\n", group_ptr->name, var_ptr2->var_array);
                  fortprintf(fd, "      end if\n");
			   } else {
				  fortprintf(fd, "      %s %% %s %% array = %s\n", group_ptr->name, var_ptr2->var_array, var_ptr2->default_value ); /* initialize field */
			   }
			} else {
               fortprintf(fd, "      %s %% %s %% isPersistent = .false.\n", group_ptr->name, var_ptr2->var_array);
               fortprintf(fd, "      %s %% %s %% isActive = .false.\n", group_ptr->name, var_ptr2->var_array);
			}

            fortprintf(fd, "      %s %% %s %% dimSizes(1) = %s %% num_%s\n", group_ptr->name, var_ptr2->var_array, group_ptr->name, var_ptr2->var_array);
            fortprintf(fd, "      %s %% %s %% dimNames(1) = \'num_%s\'\n", group_ptr->name, var_ptr2->var_array, var_ptr2->var_array);
            dimlist_ptr = var_ptr2->dimlist;
            i = 2;
            while (dimlist_ptr) {
               if (!strncmp(dimlist_ptr->dim->name_in_file, "nCells", 1024) ||
                   !strncmp(dimlist_ptr->dim->name_in_file, "nEdges", 1024) ||
                   !strncmp(dimlist_ptr->dim->name_in_file, "nVertices", 1024))
                  if (!dimlist_ptr->dim->namelist_defined) {
					 if (var_ptr2->persistence == PERSISTENT || var_ptr2->persistence == PACKAGE){
                        fortprintf(fd, "      %s %% %s %% dimSizes(%i) = %s\n", group_ptr->name, var_ptr2->var_array, i, dimlist_ptr->dim->name_in_code);
                        fortprintf(fd, "      %s %% %s %% dimNames(%i) = \'%s\'\n", group_ptr->name, var_ptr2->var_array, i, dimlist_ptr->dim->name_in_file);
					 } 
					 else {
                        fortprintf(fd, "      %s %% %s %% dimSizes(%i) = %s+1\n", group_ptr->name, var_ptr2->var_array, i, dimlist_ptr->dim->name_in_code);
                        fortprintf(fd, "      %s %% %s %% dimNames(%i) = \'%s\'\n", group_ptr->name, var_ptr2->var_array, i, dimlist_ptr->dim->name_in_file);
					 }
                  }
                  else {
                     fortprintf(fd, "      %s %% %s %% dimSizes(%i) = %s\n", group_ptr->name, var_ptr2->var_array, i, dimlist_ptr->dim->name_in_file);
                     fortprintf(fd, "      %s %% %s %% dimNames(%i) = \'%s\'\n", group_ptr->name, var_ptr2->var_array, i, dimlist_ptr->dim->name_in_file);
                  }
               else
                  if (dimlist_ptr->dim->namelist_defined) {
                     fortprintf(fd, "      %s %% %s %% dimSizes(%i) = %s\n", group_ptr->name, var_ptr2->var_array, i, dimlist_ptr->dim->name_in_file);
                     fortprintf(fd, "      %s %% %s %% dimNames(%i) = \'%s\'\n", group_ptr->name, var_ptr2->var_array, i, dimlist_ptr->dim->name_in_file);
                  }
                  else {
                     fortprintf(fd, "      %s %% %s %% dimSizes(%i) = %s\n", group_ptr->name, var_ptr2->var_array, i, dimlist_ptr->dim->name_in_code);
                     fortprintf(fd, "      %s %% %s %% dimNames(%i) = \'%s\'\n", group_ptr->name, var_ptr2->var_array, i, dimlist_ptr->dim->name_in_file);
                  }
               i++;
               dimlist_ptr = dimlist_ptr->next;
            }
            if (var_ptr2->timedim) fortprintf(fd, "      %s %% %s %% hasTimeDimension = .true.\n", group_ptr->name, var_ptr2->var_array);
            else fortprintf(fd, "      %s %% %s %% hasTimeDimension = .false.\n", group_ptr->name, var_ptr2->var_array);
            fortprintf(fd, "      nullify(%s %% %s %% prev)\n", group_ptr->name, var_ptr2->var_array);
            fortprintf(fd, "      nullify(%s %% %s %% next)\n", group_ptr->name, var_ptr2->var_array);
            fortprintf(fd, "      nullify(%s %% %s %% sendList)\n", group_ptr->name, var_ptr2->var_array);
            fortprintf(fd, "      nullify(%s %% %s %% recvList)\n", group_ptr->name, var_ptr2->var_array);
            fortprintf(fd, "      nullify(%s %% %s %% copyList)\n", group_ptr->name, var_ptr2->var_array);

            if (var_ptr2->iostreams & INPUT0) 
               fortprintf(fd, "      %s %% %s %% ioinfo %% input = .true.\n", group_ptr->name, var_ptr2->var_array);
            else
               fortprintf(fd, "      %s %% %s %% ioinfo %% input = .false.\n", group_ptr->name, var_ptr2->var_array);

            if (var_ptr2->iostreams & SFC0) 
               fortprintf(fd, "      %s %% %s %% ioinfo %% sfc = .true.\n", group_ptr->name, var_ptr2->var_array);
            else
               fortprintf(fd, "      %s %% %s %% ioinfo %% sfc = .false.\n", group_ptr->name, var_ptr2->var_array);

            if (var_ptr2->iostreams & RESTART0) 
               fortprintf(fd, "      %s %% %s %% ioinfo %% restart = .true.\n", group_ptr->name, var_ptr2->var_array);
            else
               fortprintf(fd, "      %s %% %s %% ioinfo %% restart = .false.\n", group_ptr->name, var_ptr2->var_array);

            if (var_ptr2->iostreams & OUTPUT0) 
               fortprintf(fd, "      %s %% %s %% ioinfo %% output = .true.\n", group_ptr->name, var_ptr2->var_array);
            else
               fortprintf(fd, "      %s %% %s %% ioinfo %% output = .false.\n", group_ptr->name, var_ptr2->var_array);

            fortprintf(fd, "      %s %% %s %% block => b\n", group_ptr->name, var_ptr2->var_array);
            fortprintf(fd, "\n");
         }
         else {
            fortprintf(fd, "      allocate(%s %% %s)\n", group_ptr->name, var_ptr->name_in_code);
            fortprintf(fd, "      allocate(%s %% %s %% ioinfo)\n", group_ptr->name, var_ptr->name_in_code);
            fortprintf(fd, "      %s %% %s %% fieldName = \'%s\'\n", group_ptr->name, var_ptr->name_in_code, var_ptr->name_in_file);
            fortprintf(fd, "      %s %% %s %% isVarArray = .false.\n", group_ptr->name, var_ptr->name_in_code);
            if (var_ptr->ndims > 0) {
	  		  if(var_ptr->persistence == SCRATCH){
				  fortprintf(fd, "      %s %% %s %% isPersistent = .false.\n", group_ptr->name, var_ptr->name_in_code);
				  fortprintf(fd, "      %s %% %s %% isActive = .false.\n", group_ptr->name, var_ptr->name_in_code);
				  fortprintf(fd, "      nullify(%s %% %s %% array)\n", group_ptr->name, var_ptr->name_in_code); 
			  } else if(var_ptr->persistence == PERSISTENT){
			   fortprintf(fd, "      %s %% %s %% isPersistent = .true.\n", group_ptr->name, var_ptr->name_in_code);
			   fortprintf(fd, "      %s %% %s %% isActive = .true.\n", group_ptr->name, var_ptr->name_in_code);
               fortprintf(fd, "      allocate(%s %% %s %% array(", group_ptr->name, var_ptr->name_in_code);
               dimlist_ptr = var_ptr->dimlist;
               if (!strncmp(dimlist_ptr->dim->name_in_file, "nCells", 1024) ||
                   !strncmp(dimlist_ptr->dim->name_in_file, "nEdges", 1024) ||
                   !strncmp(dimlist_ptr->dim->name_in_file, "nVertices", 1024))
                  if (!dimlist_ptr->dim->namelist_defined) fortprintf(fd, "%s + 1", dimlist_ptr->dim->name_in_code);
                  else fortprintf(fd, "%s + 1", dimlist_ptr->dim->name_in_file);
               else
                  if (dimlist_ptr->dim->namelist_defined) fortprintf(fd, "%s", dimlist_ptr->dim->name_in_file);
                  else fortprintf(fd, "%s", dimlist_ptr->dim->name_in_code);
               dimlist_ptr = dimlist_ptr->next;
               while (dimlist_ptr) {
                  if (!strncmp(dimlist_ptr->dim->name_in_file, "nCells", 1024) ||
                      !strncmp(dimlist_ptr->dim->name_in_file, "nEdges", 1024) ||
                      !strncmp(dimlist_ptr->dim->name_in_file, "nVertices", 1024))
                     if (!dimlist_ptr->dim->namelist_defined) fortprintf(fd, ", %s + 1", dimlist_ptr->dim->name_in_code);
                     else fortprintf(fd, ", %s + 1", dimlist_ptr->dim->name_in_file);
                  else
                     if (dimlist_ptr->dim->namelist_defined) fortprintf(fd, ", %s", dimlist_ptr->dim->name_in_file);
                     else fortprintf(fd, ", %s", dimlist_ptr->dim->name_in_code);
                  dimlist_ptr = dimlist_ptr->next;
               }
               fortprintf(fd, "))\n");
               fortprintf(fd, "      %s %% %s %% array = %s\n", group_ptr->name, var_ptr->name_in_code, var_ptr->default_value ); /* initialize field */
			  } else if(var_ptr->persistence == PACKAGE){
			   fortprintf(fd, "      %s %% %s %% isPersistent = .true.\n", group_ptr->name, var_ptr->name_in_code);
               fortprintf(fd, "      if (");
			   write_package_options(fd, var_ptr->package_name);
			   fortprintf(fd, ") then\n");
			   fortprintf(fd, "         %s %% %s %% isActive = .true.\n", group_ptr->name, var_ptr->name_in_code);
               fortprintf(fd, "         allocate(%s %% %s %% array(", group_ptr->name, var_ptr->name_in_code);
               dimlist_ptr = var_ptr->dimlist;
               if (!strncmp(dimlist_ptr->dim->name_in_file, "nCells", 1024) ||
                   !strncmp(dimlist_ptr->dim->name_in_file, "nEdges", 1024) ||
                   !strncmp(dimlist_ptr->dim->name_in_file, "nVertices", 1024))
                  if (!dimlist_ptr->dim->namelist_defined) fortprintf(fd, "%s + 1", dimlist_ptr->dim->name_in_code);
                  else fortprintf(fd, "%s + 1", dimlist_ptr->dim->name_in_file);
               else
                  if (dimlist_ptr->dim->namelist_defined) fortprintf(fd, "%s", dimlist_ptr->dim->name_in_file);
                  else fortprintf(fd, "%s", dimlist_ptr->dim->name_in_code);
               dimlist_ptr = dimlist_ptr->next;
               while (dimlist_ptr) {
                  if (!strncmp(dimlist_ptr->dim->name_in_file, "nCells", 1024) ||
                      !strncmp(dimlist_ptr->dim->name_in_file, "nEdges", 1024) ||
                      !strncmp(dimlist_ptr->dim->name_in_file, "nVertices", 1024))
                     if (!dimlist_ptr->dim->namelist_defined) fortprintf(fd, ", %s + 1", dimlist_ptr->dim->name_in_code);
                     else fortprintf(fd, ", %s + 1", dimlist_ptr->dim->name_in_file);
                  else
                     if (dimlist_ptr->dim->namelist_defined) fortprintf(fd, ", %s", dimlist_ptr->dim->name_in_file);
                     else fortprintf(fd, ", %s", dimlist_ptr->dim->name_in_code);
                  dimlist_ptr = dimlist_ptr->next;
               }
               fortprintf(fd, "))\n");
               if (var_ptr->vtype == INTEGER)
                  fortprintf(fd, "         %s %% %s %% array = 0\n", group_ptr->name, var_ptr->name_in_code ); /* initialize field to zero */
               else if (var_ptr->vtype == REAL)
                  fortprintf(fd, "         %s %% %s %% array = 0.0\n", group_ptr->name, var_ptr->name_in_code ); /* initialize field to zero */
               else if (var_ptr->vtype == CHARACTER)
                  fortprintf(fd, "         %s %% %s %% array = \'\'\n", group_ptr->name, var_ptr->name_in_code ); /* initialize field to zero */
			   fortprintf(fd, "      else\n");
			   fortprintf(fd, "         %s %% %s %% isActive = .false.\n", group_ptr->name, var_ptr->name_in_code);
               fortprintf(fd, "         nullify(%s %% %s %% array)\n", group_ptr->name, var_ptr->name_in_code);
			   fortprintf(fd, "      endif\n");

			  }
               dimlist_ptr = var_ptr->dimlist;
               i = 1;
               while (dimlist_ptr) {
                  if (!strncmp(dimlist_ptr->dim->name_in_file, "nCells", 1024) ||
                      !strncmp(dimlist_ptr->dim->name_in_file, "nEdges", 1024) ||
                      !strncmp(dimlist_ptr->dim->name_in_file, "nVertices", 1024))
                     if (!dimlist_ptr->dim->namelist_defined) {
						if(var_ptr->persistence == PERSISTENT || var_ptr->persistence == PACKAGE){
                          fortprintf(fd, "      %s %% %s %% dimSizes(%i) = %s\n", group_ptr->name, var_ptr->name_in_code, i, dimlist_ptr->dim->name_in_code); 
                          fortprintf(fd, "      %s %% %s %% dimNames(%i) = \'%s\'\n", group_ptr->name, var_ptr->name_in_code, i, dimlist_ptr->dim->name_in_file); 
						}
						else {
                          fortprintf(fd, "      %s %% %s %% dimSizes(%i) = %s+1\n", group_ptr->name, var_ptr->name_in_code, i, dimlist_ptr->dim->name_in_code); 
                          fortprintf(fd, "      %s %% %s %% dimNames(%i) = \'%s\'\n", group_ptr->name, var_ptr->name_in_code, i, dimlist_ptr->dim->name_in_file); 
						}
                     }
                     else {
                        fortprintf(fd, "      %s %% %s %% dimSizes(%i) = %s\n", group_ptr->name, var_ptr->name_in_code, i, dimlist_ptr->dim->name_in_file); 
                        fortprintf(fd, "      %s %% %s %% dimNames(%i) = \'%s\'\n", group_ptr->name, var_ptr->name_in_code, i, dimlist_ptr->dim->name_in_file); 
                     }
                  else
                     if (dimlist_ptr->dim->namelist_defined) {
                        fortprintf(fd, "      %s %% %s %% dimSizes(%i) = %s\n", group_ptr->name, var_ptr->name_in_code, i, dimlist_ptr->dim->name_in_file); 
                        fortprintf(fd, "      %s %% %s %% dimNames(%i) = \'%s\'\n", group_ptr->name, var_ptr->name_in_code, i, dimlist_ptr->dim->name_in_file); 
                     }
                     else {
                        fortprintf(fd, "      %s %% %s %% dimSizes(%i) = %s\n", group_ptr->name, var_ptr->name_in_code, i, dimlist_ptr->dim->name_in_code); 
                        fortprintf(fd, "      %s %% %s %% dimNames(%i) = \'%s\'\n", group_ptr->name, var_ptr->name_in_code, i, dimlist_ptr->dim->name_in_file); 
                     }
                  i++;
                  dimlist_ptr = dimlist_ptr->next;
               }
			}

            if (var_ptr->timedim) fortprintf(fd, "      %s %% %s %% hasTimeDimension = .true.\n", group_ptr->name, var_ptr->name_in_code);
            else fortprintf(fd, "      %s %% %s %% hasTimeDimension = .false.\n", group_ptr->name, var_ptr->name_in_code);
            fortprintf(fd, "      nullify(%s %% %s %% prev)\n", group_ptr->name, var_ptr->name_in_code);
            fortprintf(fd, "      nullify(%s %% %s %% next)\n", group_ptr->name, var_ptr->name_in_code);
            fortprintf(fd, "      nullify(%s %% %s %% sendList)\n", group_ptr->name, var_ptr->name_in_code);
            fortprintf(fd, "      nullify(%s %% %s %% recvList)\n", group_ptr->name, var_ptr->name_in_code);
            fortprintf(fd, "      nullify(%s %% %s %% copyList)\n", group_ptr->name, var_ptr->name_in_code);

            if (var_ptr->iostreams & INPUT0) 
               fortprintf(fd, "      %s %% %s %% ioinfo %% input = .true.\n", group_ptr->name, var_ptr->name_in_code);
            else
               fortprintf(fd, "      %s %% %s %% ioinfo %% input = .false.\n", group_ptr->name, var_ptr->name_in_code);

            if (var_ptr->iostreams & SFC0) 
               fortprintf(fd, "      %s %% %s %% ioinfo %% sfc = .true.\n", group_ptr->name, var_ptr->name_in_code);
            else
               fortprintf(fd, "      %s %% %s %% ioinfo %% sfc = .false.\n", group_ptr->name, var_ptr->name_in_code);

            if (var_ptr->iostreams & RESTART0) 
               fortprintf(fd, "      %s %% %s %% ioinfo %% restart = .true.\n", group_ptr->name, var_ptr->name_in_code);
            else
               fortprintf(fd, "      %s %% %s %% ioinfo %% restart = .false.\n", group_ptr->name, var_ptr->name_in_code);

            if (var_ptr->iostreams & OUTPUT0) 
               fortprintf(fd, "      %s %% %s %% ioinfo %% output = .true.\n", group_ptr->name, var_ptr->name_in_code);
            else
               fortprintf(fd, "      %s %% %s %% ioinfo %% output = .false.\n", group_ptr->name, var_ptr->name_in_code);

            fortprintf(fd, "      %s %% %s %% block => b\n", group_ptr->name, var_ptr->name_in_code);
            fortprintf(fd, "\n");

            var_list_ptr = var_list_ptr->next;
         }
      }

      fortprintf(fd, "   end subroutine mpas_allocate_%s\n\n\n", group_ptr->name);
      group_ptr = group_ptr->next;
   }
   fclose(fd);
   
   /* Definitions of deallocate subroutines */
   fd = fopen("group_dealloc_routines.inc", "w");
   group_ptr = groups;
   while (group_ptr) {
      fortprintf(fd, "   subroutine mpas_deallocate_%s(%s)\n", group_ptr->name, group_ptr->name);
      fortprintf(fd, "\n");
      fortprintf(fd, "      implicit none\n");
      fortprintf(fd, "\n");
      fortprintf(fd, "      type (%s_type), intent(inout) :: %s\n", group_ptr->name, group_ptr->name);
      fortprintf(fd, "\n");

      var_list_ptr = group_ptr->vlist;
      while (var_list_ptr) {
         var_ptr = var_list_ptr->var;
         if (strncmp(var_ptr->var_array, "-", 1024) != 0) {
            memcpy(var_array, var_ptr->var_array, 1024);
            memcpy(array_class, var_ptr->array_class, 1024);
            i = 0;
            while (var_list_ptr && strncmp(var_array, var_list_ptr->var->var_array, 1024) == 0) {
               i++;
               var_list_ptr2 = var_list_ptr;
               var_list_ptr = var_list_ptr->next;
            }
            fortprintf(fd, "      if(associated(%s %% %s %% array)) then\n", group_ptr->name, var_list_ptr2->var->var_array);
            fortprintf(fd, "         deallocate(%s %% %s %% array)\n", group_ptr->name, var_list_ptr2->var->var_array);
            fortprintf(fd, "      end if\n");
            fortprintf(fd, "      deallocate(%s %% %s %% ioinfo)\n", group_ptr->name, var_list_ptr2->var->var_array);
            fortprintf(fd, "      call mpas_deallocate_attlist(%s %% %s %% attList)\n", group_ptr->name, var_list_ptr2->var->var_array);
            fortprintf(fd, "      deallocate(%s %% %s)\n\n", group_ptr->name, var_list_ptr2->var->var_array);
         }
         else {
            if (var_ptr->ndims > 0) {
               fortprintf(fd, "      if(associated(%s %% %s %% array)) then\n", group_ptr->name, var_ptr->name_in_code);
               fortprintf(fd, "         deallocate(%s %% %s %% array)\n", group_ptr->name, var_ptr->name_in_code);
               fortprintf(fd, "      end if\n");
               fortprintf(fd, "      deallocate(%s %% %s %% ioinfo)\n", group_ptr->name, var_ptr->name_in_code);
               fortprintf(fd, "      call mpas_deallocate_attlist(%s %% %s %% attList)\n", group_ptr->name, var_ptr->name_in_code);
               fortprintf(fd, "      deallocate(%s %% %s)\n\n", group_ptr->name, var_ptr->name_in_code);
            }
            else {
               fortprintf(fd, "      deallocate(%s %% %s %% ioinfo)\n", group_ptr->name, var_ptr->name_in_code);
               fortprintf(fd, "      call mpas_deallocate_attlist(%s %% %s %% attList)\n", group_ptr->name, var_ptr->name_in_code);
               fortprintf(fd, "      deallocate(%s %% %s)\n\n", group_ptr->name, var_ptr->name_in_code);
            }
            var_list_ptr = var_list_ptr->next;
         }
      }

      fortprintf(fd, "   end subroutine mpas_deallocate_%s\n\n\n", group_ptr->name);
      group_ptr = group_ptr->next;
   }
   fclose(fd);

   /* Definitions of copy subroutines */
   fd = fopen("group_copy_routines.inc", "w");
   group_ptr = groups;
   while (group_ptr) {
      fortprintf(fd, "   subroutine mpas_copy_%s(dest, src)\n", group_ptr->name);
      fortprintf(fd, "\n");
      fortprintf(fd, "      implicit none\n");
      fortprintf(fd, "\n");
      fortprintf(fd, "      type (%s_type), intent(in) :: src\n", group_ptr->name);
      fortprintf(fd, "      type (%s_type), intent(inout) :: dest\n", group_ptr->name);
      fortprintf(fd, "\n");
      var_list_ptr = group_ptr->vlist;
      while (var_list_ptr) {
         var_ptr = var_list_ptr->var;
         if (strncmp(var_ptr->var_array, "-", 1024) != 0) {
            memcpy(var_array, var_ptr->var_array, 1024);
            memcpy(array_class, var_ptr->array_class, 1024);
            i = 0;
            while (var_list_ptr && strncmp(var_array, var_list_ptr->var->var_array, 1024) == 0) {
               i++;
               var_list_ptr2 = var_list_ptr;
               var_list_ptr = var_list_ptr->next;
            }
            var_ptr2 = var_list_ptr2->var;
            if (var_ptr2->ndims > 0)  {
			   fortprintf(fd, "      if (associated(dest %% %s %% array) .and. associated(src %% %s %% array)) then\n", var_ptr2->var_array, var_ptr2->var_array);
               fortprintf(fd, "         dest %% %s %% array = src %% %s %% array\n", var_ptr2->var_array, var_ptr2->var_array);
			   fortprintf(fd, "      end if\n");
            } else {
               fortprintf(fd, "      dest %% %s %% scalar = src %% %s %% scalar\n", var_ptr2->var_array, var_ptr2->var_array);
			}
         }
         else {
            if (var_ptr->ndims > 0) {
			   fortprintf(fd, "      if (associated(dest %% %s %% array) .and. associated(src %% %s %% array)) then\n", var_ptr->name_in_code, var_ptr->name_in_code);
               fortprintf(fd, "         dest %% %s %% array = src %% %s %% array\n", var_ptr->name_in_code, var_ptr->name_in_code);
			   fortprintf(fd, "      end if\n");
            } else {
               fortprintf(fd, "      dest %% %s %% scalar = src %% %s %% scalar\n", var_ptr->name_in_code, var_ptr->name_in_code);
			}
            var_list_ptr = var_list_ptr->next;
         }
      }
      fortprintf(fd, "\n");
      fortprintf(fd, "   end subroutine mpas_copy_%s\n\n\n", group_ptr->name);
      group_ptr = group_ptr->next;
   }
   fclose(fd);

   /* Definitions of shift_time_level subroutines */
   fd = fopen("group_shift_level_routines.inc", "w");
   group_ptr = groups;
   while (group_ptr) {
      if (group_ptr->vlist != NULL && group_ptr->ntime_levs > 1) {
         fortprintf(fd, "   subroutine mpas_shift_time_levels_%s(%s)\n", group_ptr->name, group_ptr->name);
         fortprintf(fd, "\n");
         fortprintf(fd, "      implicit none\n");
         fortprintf(fd, "\n");
         fortprintf(fd, "      type (%s_multilevel_type), intent(inout) :: %s\n", group_ptr->name, group_ptr->name);
         fortprintf(fd, "\n");
         fortprintf(fd, "      integer :: i\n");
         fortprintf(fd, "      real (kind=RKIND) :: real0d\n");
         fortprintf(fd, "      real (kind=RKIND), dimension(:), pointer :: real1d\n");
         fortprintf(fd, "      real (kind=RKIND), dimension(:,:), pointer :: real2d\n");
         fortprintf(fd, "      real (kind=RKIND), dimension(:,:,:), pointer :: real3d\n");
         fortprintf(fd, "      integer :: int0d\n");
         fortprintf(fd, "      integer, dimension(:), pointer :: int1d\n");
         fortprintf(fd, "      integer, dimension(:,:), pointer :: int2d\n");
         fortprintf(fd, "      integer, dimension(:,:,:), pointer :: int3d\n");
         fortprintf(fd, "      character (len=64) :: char0d\n");
         fortprintf(fd, "      character (len=64), dimension(:), pointer :: char1d\n");
         fortprintf(fd, "\n");
         var_list_ptr = group_ptr->vlist;
         while (var_list_ptr) {
            var_ptr = var_list_ptr->var;

            if (strncmp(var_ptr->var_array, "-", 1024) != 0) 
            {
               if (var_ptr->vtype == INTEGER) sprintf(type_str, "int%id", var_ptr->ndims+1); 
               else if (var_ptr->vtype == REAL) sprintf(type_str, "real%id", var_ptr->ndims+1); 
               else if (var_ptr->vtype == CHARACTER) sprintf(type_str, "char%id", var_ptr->ndims+1); 

               memcpy(var_array, var_ptr->var_array, 1024);

               while (var_list_ptr && strncmp(var_array, var_list_ptr->var->var_array, 1024) == 0)
               {
                  var_list_ptr2 = var_list_ptr;
                  var_list_ptr = var_list_ptr->next;
               }
               var_ptr2 = var_list_ptr2->var;

               fortprintf(fd, "      %s => %s %% time_levs(1) %% %s %% %s %% array\n", type_str, group_ptr->name, group_ptr->name, var_ptr2->var_array);

               fortprintf(fd, "      do i=1,%s %% nTimeLevels-1\n", group_ptr->name);
               fortprintf(fd, "         %s %% time_levs(i) %% %s %% %s %% array => %s %% time_levs(i+1) %% %s %% %s %% array\n", group_ptr->name, group_ptr->name, var_ptr2->var_array, group_ptr->name, group_ptr->name, var_ptr2->var_array);
               fortprintf(fd, "      end do\n");

               fortprintf(fd, "      %s %% time_levs(%s %% nTimeLevels) %% %s %% %s %% array=> %s\n\n", group_ptr->name, group_ptr->name, group_ptr->name, var_ptr2->var_array, type_str);
            }
            else {

               if (var_ptr->vtype == INTEGER) sprintf(type_str, "int%id", var_ptr->ndims); 
               else if (var_ptr->vtype == REAL) sprintf(type_str, "real%id", var_ptr->ndims); 
               else if (var_ptr->vtype == CHARACTER) sprintf(type_str, "char%id", var_ptr->ndims); 

               if (var_ptr->ndims > 0) 
                  fortprintf(fd, "      %s => %s %% time_levs(1) %% %s %% %s %% array\n", type_str, group_ptr->name, group_ptr->name, var_ptr->name_in_code);
               else
                  fortprintf(fd, "      %s = %s %% time_levs(1) %% %s %% %s %% scalar\n", type_str, group_ptr->name, group_ptr->name, var_ptr->name_in_code);

               fortprintf(fd, "      do i=1,%s %% nTimeLevels-1\n", group_ptr->name);
               if (var_ptr->ndims > 0) 
                  fortprintf(fd, "         %s %% time_levs(i) %% %s %% %s %% array => %s %% time_levs(i+1) %% %s %% %s %% array\n", group_ptr->name, group_ptr->name, var_ptr->name_in_code, group_ptr->name, group_ptr->name, var_ptr->name_in_code);
               else
                  fortprintf(fd, "         %s %% time_levs(i) %% %s %% %s %% scalar = %s %% time_levs(i+1) %% %s %% %s %% scalar\n", group_ptr->name, group_ptr->name, var_ptr->name_in_code, group_ptr->name, group_ptr->name, var_ptr->name_in_code);
               fortprintf(fd, "      end do\n");

               if (var_ptr->ndims > 0) 
                  fortprintf(fd, "      %s %% time_levs(%s %% nTimeLevels) %% %s %% %s %% array=> %s\n\n", group_ptr->name, group_ptr->name, group_ptr->name, var_ptr->name_in_code, type_str);
               else
                  fortprintf(fd, "      %s %% time_levs(%s %% nTimeLevels) %% %s %% %s %% scalar = %s\n\n", group_ptr->name, group_ptr->name, group_ptr->name, var_ptr->name_in_code, type_str);

               var_list_ptr = var_list_ptr->next;
            }
         }
         fortprintf(fd, "\n");
         fortprintf(fd, "   end subroutine mpas_shift_time_levels_%s\n\n\n", group_ptr->name);
      }
      group_ptr = group_ptr->next;
   }
   fclose(fd);
   

   /* Definitions of deallocate subroutines */
   fd = fopen("field_links.inc", "w");

   /* subroutine to call link subroutine for every field type */
   fortprintf(fd, "   subroutine mpas_create_field_links(b)\n\n");

   fortprintf(fd, "      implicit none\n");

   fortprintf(fd, "      type (block_type), pointer :: b\n");
   fortprintf(fd, "      type (block_type), pointer :: prev, next\n\n");

   fortprintf(fd, "      if (associated(b %% prev)) then\n");
   fortprintf(fd, "         prev => b %% prev\n");
   fortprintf(fd, "      else\n");
   fortprintf(fd, "         nullify(prev)\n");
   fortprintf(fd, "      end if\n");
   fortprintf(fd, "      if (associated(b %% next)) then\n");
   fortprintf(fd, "         next => b %% next\n");
   fortprintf(fd, "      else\n");
   fortprintf(fd, "         nullify(next)\n");
   fortprintf(fd, "      end if\n\n");

   group_ptr = groups;
   while (group_ptr)
   {
      var_list_ptr = group_ptr->vlist;

      if (!group_ptr->vlist) {
         group_ptr = group_ptr->next;
         continue;
      }

      if (group_ptr->ntime_levs > 1) {
         for(i=1; i<=group_ptr->ntime_levs; i++) {
            fortprintf(fd, "      if (associated(next) .and. associated(prev)) then\n");	
            fortprintf(fd, "         call mpas_create_%s_links(b %% %s %% time_levs(%i) %% %s, ", group_ptr->name, group_ptr->name, i, group_ptr->name, i);
            fortprintf(fd, " prev = prev %% %s %% time_levs(%i) %% %s,", group_ptr->name, i, group_ptr->name);
            fortprintf(fd, " next = next %% %s %% time_levs(%i) %% %s)\n", group_ptr->name, i, group_ptr->name);
            fortprintf(fd, "      else if (associated(next)) then\n");	
            fortprintf(fd, "         call mpas_create_%s_links(b %% %s %% time_levs(%i) %% %s, next = next %% %s %% time_levs(%i) %% %s)\n", group_ptr->name, group_ptr->name, i, group_ptr->name, group_ptr->name, i, group_ptr->name);
            fortprintf(fd, "      else if (associated(prev)) then\n");	
            fortprintf(fd, "         call mpas_create_%s_links(b %% %s %% time_levs(%i) %% %s, prev = prev %% %s %% time_levs(%i) %% %s)\n", group_ptr->name, group_ptr->name, i, group_ptr->name, group_ptr->name, i, group_ptr->name);
            fortprintf(fd, "      else\n");
            fortprintf(fd, "         call mpas_create_%s_links(b %% %s %% time_levs(%i) %% %s)\n", group_ptr->name, group_ptr->name, i, group_ptr->name);
            fortprintf(fd, "      end if\n\n");
         }	
      }
      else {
         fortprintf(fd, "      if (associated(next) .and. associated(prev)) then\n");	
         fortprintf(fd, "         call mpas_create_%s_links(b %% %s, prev = prev %% %s, next = next %% %s)\n", group_ptr->name, group_ptr->name, group_ptr->name, group_ptr->name); 
         fortprintf(fd, "      else if (associated(next)) then\n");	
         fortprintf(fd, "         call mpas_create_%s_links(b %% %s, next = next %% %s)\n", group_ptr->name, group_ptr->name, group_ptr->name); 
         fortprintf(fd, "      else if (associated(prev)) then\n");	
         fortprintf(fd, "         call mpas_create_%s_links(b %% %s, prev = prev %% %s)\n", group_ptr->name, group_ptr->name, group_ptr->name); 
         fortprintf(fd, "      else\n");
         fortprintf(fd, "         call mpas_create_%s_links(b %% %s)\n", group_ptr->name, group_ptr->name); 
         fortprintf(fd, "      end if\n\n");
      }

      group_ptr = group_ptr->next;
   }
   fortprintf(fd, "   end subroutine mpas_create_field_links\n\n\n");

   /* subroutines for linking specific field type */
   group_ptr = groups;

   while (group_ptr) {
      fortprintf(fd, "      subroutine mpas_create_%s_links(%s, prev, next)\n\n", group_ptr->name, group_ptr->name); 
      fortprintf(fd, "         implicit none\n");
      fortprintf(fd, "         type (%s_type), pointer :: %s\n", group_ptr->name, group_ptr->name);
	  fortprintf(fd, "         type (%s_type), pointer, optional :: prev, next\n", group_ptr->name);

      var_list_ptr = group_ptr->vlist;
      while (var_list_ptr) {
         var_ptr = var_list_ptr->var;
         if (strncmp(var_ptr->var_array, "-", 1024) != 0) {
            memcpy(var_array, var_ptr->var_array, 1024);
            memcpy(array_class, var_ptr->array_class, 1024);
            while (var_list_ptr && strncmp(var_array, var_list_ptr->var->var_array, 1024) == 0) {
               var_list_ptr2 = var_list_ptr;
               var_list_ptr = var_list_ptr->next;
            }
            var_ptr2 = var_list_ptr2->var;
            get_outer_dim(var_ptr2, outer_dim);
            
               if (strncmp("nCells",outer_dim,1024) == 0) {
                  fortprintf(fd, "         %s %% %s %% sendList => %s %% %s %% block %% parinfo %% cellsToSend\n", group_ptr->name, var_ptr2->var_array, group_ptr->name, var_ptr2->var_array);
                  fortprintf(fd, "         %s %% %s %% recvList => %s %% %s %% block %% parinfo %% cellsToRecv\n", group_ptr->name, var_ptr2->var_array, group_ptr->name, var_ptr2->var_array);
                  fortprintf(fd, "         %s %% %s %% copyList => %s %% %s %% block %% parinfo %% cellsToCopy\n", group_ptr->name, var_ptr2->var_array, group_ptr->name, var_ptr2->var_array);
				  fortprintf(fd, "         if(present(prev)) then\n");
				  fortprintf(fd, "           %s %% %s %% prev => prev %% %s\n", group_ptr->name, var_ptr2->var_array, var_ptr2->var_array);
				  fortprintf(fd, "         else\n");
				  fortprintf(fd, "           nullify(%s %% %s %% prev)\n", group_ptr->name, var_ptr2->var_array);
				  fortprintf(fd, "         end if\n");
				  fortprintf(fd, "         if(present(next)) then\n");
				  fortprintf(fd, "           %s %% %s %% next => next %% %s\n", group_ptr->name, var_ptr2->var_array, var_ptr2->var_array);
				  fortprintf(fd, "         else\n");
				  fortprintf(fd, "           nullify(%s %% %s %% next)\n", group_ptr->name, var_ptr2->var_array);
				  fortprintf(fd, "         end if\n\n");
               }
               else if (strncmp("nEdges",outer_dim,1024) == 0) {
                  fortprintf(fd, "         %s %% %s %% sendList => %s %% %s %% block %% parinfo %% edgesToSend\n", group_ptr->name, var_ptr2->var_array, group_ptr->name, var_ptr2->var_array);
                  fortprintf(fd, "         %s %% %s %% recvList => %s %% %s %% block %% parinfo %% edgesToRecv\n", group_ptr->name, var_ptr2->var_array, group_ptr->name, var_ptr2->var_array);
                  fortprintf(fd, "         %s %% %s %% copyList => %s %% %s %% block %% parinfo %% edgesToCopy\n", group_ptr->name, var_ptr2->var_array, group_ptr->name, var_ptr2->var_array);
				  fortprintf(fd, "         if(present(prev)) then\n");
				  fortprintf(fd, "           %s %% %s %% prev => prev %% %s\n", group_ptr->name, var_ptr2->var_array, var_ptr2->var_array);
				  fortprintf(fd, "         else\n");
				  fortprintf(fd, "           nullify(%s %% %s %% prev)\n", group_ptr->name, var_ptr2->var_array);
				  fortprintf(fd, "         end if\n");
				  fortprintf(fd, "         if(present(next)) then\n");
				  fortprintf(fd, "           %s %% %s %% next => next %% %s\n", group_ptr->name, var_ptr2->var_array, var_ptr2->var_array);
				  fortprintf(fd, "         else\n");
				  fortprintf(fd, "           nullify(%s %% %s %% next)\n", group_ptr->name, var_ptr2->var_array);
				  fortprintf(fd, "         end if\n\n");
               }
               else if (strncmp("nVertices",outer_dim,1024) == 0) {
                  fortprintf(fd, "         %s %% %s %% sendList => %s %% %s %% block %% parinfo %% verticesToSend\n", group_ptr->name, var_ptr2->var_array, group_ptr->name, var_ptr2->var_array);
                  fortprintf(fd, "         %s %% %s %% recvList => %s %% %s %% block %% parinfo %% verticesToRecv\n", group_ptr->name, var_ptr2->var_array, group_ptr->name, var_ptr2->var_array);
                  fortprintf(fd, "         %s %% %s %% copyList => %s %% %s %% block %% parinfo %% verticesToCopy\n", group_ptr->name, var_ptr2->var_array, group_ptr->name, var_ptr2->var_array);
				  fortprintf(fd, "         if(present(prev)) then\n");
				  fortprintf(fd, "           %s %% %s %% prev => prev %% %s\n", group_ptr->name, var_ptr2->var_array, var_ptr2->var_array);
				  fortprintf(fd, "         else\n");
				  fortprintf(fd, "           nullify(%s %% %s %% prev)\n", group_ptr->name, var_ptr2->var_array);
				  fortprintf(fd, "         end if\n");
				  fortprintf(fd, "         if(present(next)) then\n");
				  fortprintf(fd, "           %s %% %s %% next => next %% %s\n", group_ptr->name, var_ptr2->var_array, var_ptr2->var_array);
				  fortprintf(fd, "         else\n");
				  fortprintf(fd, "           nullify(%s %% %s %% next)\n", group_ptr->name, var_ptr2->var_array);
				  fortprintf(fd, "         end if\n\n");
               } else {
				  fortprintf(fd, "         nullify(%s %% %s %% sendList)\n", group_ptr->name, var_ptr2->var_array);
				  fortprintf(fd, "         nullify(%s %% %s %% recvList)\n", group_ptr->name, var_ptr2->var_array);
				  fortprintf(fd, "         nullify(%s %% %s %% copyList)\n", group_ptr->name, var_ptr2->var_array);
				  fortprintf(fd, "         if(present(prev)) then\n");
				  fortprintf(fd, "           %s %% %s %% prev => prev %% %s\n", group_ptr->name, var_ptr2->var_array, var_ptr2->var_array);
				  fortprintf(fd, "         else\n");
				  fortprintf(fd, "           nullify(%s %% %s %% prev)\n", group_ptr->name, var_ptr2->var_array);
				  fortprintf(fd, "         end if\n");
				  fortprintf(fd, "         if(present(next)) then\n");
				  fortprintf(fd, "           %s %% %s %% next => next %% %s\n", group_ptr->name, var_ptr2->var_array, var_ptr2->var_array);
				  fortprintf(fd, "         else\n");
				  fortprintf(fd, "           nullify(%s %% %s %% next)\n", group_ptr->name, var_ptr2->var_array);
				  fortprintf(fd, "         end if\n\n");

			   }
            fortprintf(fd, "\n");
         }
         else 
         {
	    if (var_ptr->ndims > 0)
	    {
               get_outer_dim(var_ptr, outer_dim);
               
               if (strncmp("nCells",outer_dim,1024) == 0) {
                  fortprintf(fd, "         %s %% %s %% sendList => %s %% %s %% block %% parinfo %% cellsToSend\n", group_ptr->name, var_ptr->name_in_code, group_ptr->name, var_ptr->name_in_code);
                  fortprintf(fd, "         %s %% %s %% recvList => %s %% %s %% block %% parinfo %% cellsToRecv\n", group_ptr->name, var_ptr->name_in_code, group_ptr->name, var_ptr->name_in_code);
                  fortprintf(fd, "         %s %% %s %% copyList => %s %% %s %% block %% parinfo %% cellsToCopy\n", group_ptr->name, var_ptr->name_in_code, group_ptr->name, var_ptr->name_in_code);
				  fortprintf(fd, "         if(present(prev)) then\n");
				  fortprintf(fd, "           %s %% %s %% prev => prev %% %s\n", group_ptr->name, var_ptr->name_in_code, var_ptr->name_in_code);
				  fortprintf(fd, "         else\n");
				  fortprintf(fd, "           nullify(%s %% %s %% prev)\n", group_ptr->name, var_ptr->name_in_code);
				  fortprintf(fd, "         end if\n");
				  fortprintf(fd, "         if(present(next)) then\n");
				  fortprintf(fd, "           %s %% %s %% next => next %% %s\n", group_ptr->name, var_ptr->name_in_code, var_ptr->name_in_code);
				  fortprintf(fd, "         else\n");
				  fortprintf(fd, "           nullify(%s %% %s %% next)\n", group_ptr->name, var_ptr->name_in_code);
				  fortprintf(fd, "         end if\n\n");
               }
               else if (strncmp("nEdges",outer_dim,1024) == 0) {
                  fortprintf(fd, "         %s %% %s %% sendList => %s %% %s %% block %% parinfo %% edgesToSend\n", group_ptr->name, var_ptr->name_in_code, group_ptr->name, var_ptr->name_in_code);
                  fortprintf(fd, "         %s %% %s %% recvList => %s %% %s %% block %% parinfo %% edgesToRecv\n", group_ptr->name, var_ptr->name_in_code, group_ptr->name, var_ptr->name_in_code);
                  fortprintf(fd, "         %s %% %s %% copyList => %s %% %s %% block %% parinfo %% edgesToCopy\n", group_ptr->name, var_ptr->name_in_code, group_ptr->name, var_ptr->name_in_code);
				  fortprintf(fd, "         if(present(prev)) then\n");
				  fortprintf(fd, "           %s %% %s %% prev => prev %% %s\n", group_ptr->name, var_ptr->name_in_code, var_ptr->name_in_code);
				  fortprintf(fd, "         else\n");
				  fortprintf(fd, "           nullify(%s %% %s %% prev)\n", group_ptr->name, var_ptr->name_in_code);
				  fortprintf(fd, "         end if\n");
				  fortprintf(fd, "         if(present(next)) then\n");
				  fortprintf(fd, "           %s %% %s %% next => next %% %s\n", group_ptr->name, var_ptr->name_in_code, var_ptr->name_in_code);
				  fortprintf(fd, "         else\n");
				  fortprintf(fd, "           nullify(%s %% %s %% next)\n", group_ptr->name, var_ptr->name_in_code);
				  fortprintf(fd, "         end if\n\n");
               }
               else if (strncmp("nVertices",outer_dim,1024) == 0) {
                  fortprintf(fd, "         %s %% %s %% sendList => %s %% %s %% block %% parinfo %% verticesToSend\n", group_ptr->name, var_ptr->name_in_code, group_ptr->name, var_ptr->name_in_code);
                  fortprintf(fd, "         %s %% %s %% recvList => %s %% %s %% block %% parinfo %% verticesToRecv\n", group_ptr->name, var_ptr->name_in_code, group_ptr->name, var_ptr->name_in_code);
                  fortprintf(fd, "         %s %% %s %% copyList => %s %% %s %% block %% parinfo %% verticesToCopy\n", group_ptr->name, var_ptr->name_in_code, group_ptr->name, var_ptr->name_in_code);
				  fortprintf(fd, "         if(present(prev)) then\n");
				  fortprintf(fd, "           %s %% %s %% prev => prev %% %s\n", group_ptr->name, var_ptr->name_in_code, var_ptr->name_in_code);
				  fortprintf(fd, "         else\n");
				  fortprintf(fd, "           nullify(%s %% %s %% prev)\n", group_ptr->name, var_ptr->name_in_code);
				  fortprintf(fd, "         end if\n");
				  fortprintf(fd, "         if(present(next)) then\n");
				  fortprintf(fd, "           %s %% %s %% next => next %% %s\n", group_ptr->name, var_ptr->name_in_code, var_ptr->name_in_code);
				  fortprintf(fd, "         else\n");
				  fortprintf(fd, "           nullify(%s %% %s %% next)\n", group_ptr->name, var_ptr->name_in_code);
				  fortprintf(fd, "         end if\n\n");
               } else {
                  fortprintf(fd, "         nullify(%s %% %s %% sendList)\n", group_ptr->name, var_ptr->name_in_code);
                  fortprintf(fd, "         nullify(%s %% %s %% recvList)\n", group_ptr->name, var_ptr->name_in_code);
                  fortprintf(fd, "         nullify(%s %% %s %% copyList)\n", group_ptr->name, var_ptr->name_in_code);
				  fortprintf(fd, "         if(present(prev)) then\n");
				  fortprintf(fd, "           %s %% %s %% prev => prev %% %s\n", group_ptr->name, var_ptr->name_in_code, var_ptr->name_in_code);
				  fortprintf(fd, "         else\n");
				  fortprintf(fd, "           nullify(%s %% %s %% prev)\n", group_ptr->name, var_ptr->name_in_code);
				  fortprintf(fd, "         end if\n");
				  fortprintf(fd, "         if(present(next)) then\n");
				  fortprintf(fd, "           %s %% %s %% next => next %% %s\n", group_ptr->name, var_ptr->name_in_code, var_ptr->name_in_code);
				  fortprintf(fd, "         else\n");
				  fortprintf(fd, "           nullify(%s %% %s %% next)\n", group_ptr->name, var_ptr->name_in_code);
				  fortprintf(fd, "         end if\n\n");
			   }
               fortprintf(fd, "\n");
	    }
            var_list_ptr = var_list_ptr->next;
	 }
      }
     
      fortprintf(fd, "      end subroutine mpas_create_%s_links\n\n\n", group_ptr->name); 

      group_ptr = group_ptr->next;
   }
   fclose(fd);
 
}


void gen_reads(struct group_list * groups, struct variable * vars, struct dimension * dims)
{
   struct variable * var_ptr;
   struct variable_list * var_list_ptr, *var_list_ptr2;
   struct dimension * dim_ptr;
   struct dimension_list * dimlist_ptr, * lastdim;
   struct group_list * group_ptr;
   struct dtable * dictionary;
   FILE * fd, *fd2;
   char vtype[5];
   char fname[32];
   char var_array[1024];
   char struct_deref[1024];
   char * cp1, * cp2;
   int i, j;
   int ivtype;


   fd = fopen("add_input_fields.inc", "w");

   group_ptr = groups;
   while (group_ptr) {
      var_list_ptr = group_ptr->vlist;
      while (var_list_ptr) {
         var_ptr = var_list_ptr->var;

         if (group_ptr->ntime_levs > 1)
            snprintf(struct_deref, 1024, "blocklist %% %s %% time_levs(1) %% %s", group_ptr->name, group_ptr->name);
         else
            snprintf(struct_deref, 1024, "blocklist %% %s", group_ptr->name);
         
         if (strncmp(var_ptr->var_array, "-", 1024) != 0) {
            fortprintf(fd, "      if ((%s %% %s %% ioinfo %% input .and. input_obj %% stream == STREAM_INPUT) .or. &\n", struct_deref, var_ptr->var_array);
            fortprintf(fd, "          (%s %% %s %% ioinfo %% restart .and. input_obj %% stream == STREAM_RESTART) .or. &\n", struct_deref, var_ptr->var_array);
            fortprintf(fd, "          (%s %% %s %% ioinfo %% sfc .and. input_obj %% stream == STREAM_SFC)) then\n", struct_deref, var_ptr->var_array);
            memcpy(var_array, var_ptr->var_array, 1024);
/*            fortprintf(fd, "         write(0,*) \'adding input field %s\'\n", var_ptr->var_array); */
            fortprintf(fd, "         call MPAS_streamAddField(input_obj %% io_stream, %s %% %s, nferr)\n", struct_deref, var_ptr->var_array);
            while (var_list_ptr && strncmp(var_array, var_list_ptr->var->var_array, 1024) == 0) {
			   var_list_ptr2 = var_list_ptr;
               var_list_ptr = var_list_ptr->next;
            }
			var_list_ptr = var_list_ptr2;
         }
         else {
            fortprintf(fd, "      if ((%s %% %s %% ioinfo %% input .and. input_obj %% stream == STREAM_INPUT) .or. &\n", struct_deref, var_ptr->name_in_code);
            fortprintf(fd, "          (%s %% %s %% ioinfo %% restart .and. input_obj %% stream == STREAM_RESTART) .or. &\n", struct_deref, var_ptr->name_in_code);
            fortprintf(fd, "          (%s %% %s %% ioinfo %% sfc .and. input_obj %% stream == STREAM_SFC)) then\n", struct_deref, var_ptr->name_in_code);
/*            fortprintf(fd, "         write(0,*) \'adding input field %s\'\n", var_ptr->name_in_code); */
            fortprintf(fd, "         call MPAS_streamAddField(input_obj %% io_stream, %s %% %s, nferr)\n", struct_deref, var_ptr->name_in_code);
         }
   
         fortprintf(fd, "      end if\n\n");

         if (var_list_ptr) var_list_ptr = var_list_ptr->next;
      }
      group_ptr = group_ptr->next;
   }

   fclose(fd);


   fd = fopen("exchange_input_field_halos.inc", "w");
   fd2 = fopen("non_decomp_copy_input_fields.inc", "w");

   group_ptr = groups;
   while (group_ptr) {
      var_list_ptr = group_ptr->vlist;
      while (var_list_ptr) {
         var_ptr = var_list_ptr->var;

         dimlist_ptr = var_ptr->dimlist;
         i = 1;
		 if(var_ptr->persistence == PERSISTENT || var_ptr->persistence == PACKAGE){
         while (dimlist_ptr) {
            if (i == var_ptr->ndims) { 

                  if (group_ptr->ntime_levs > 1) {
                     snprintf(struct_deref, 1024, "domain %% blocklist %% %s %% time_levs(1) %% %s", group_ptr->name, group_ptr->name);
				  } else {
                     snprintf(struct_deref, 1024, "domain %% blocklist %% %s", group_ptr->name);
				  }

               if (!strncmp(dimlist_ptr->dim->name_in_file, "nCells", 1024) ||
                   !strncmp(dimlist_ptr->dim->name_in_file, "nEdges", 1024) ||
                   !strncmp(dimlist_ptr->dim->name_in_file, "nVertices", 1024)) {
                  
                  if (strncmp(var_ptr->var_array, "-", 1024) != 0) {
                     fortprintf(fd, "      if (%s %% %s %% isPersistent .and. %s %% %s %% isActive) then\n", struct_deref, var_ptr->var_array, struct_deref, var_ptr->var_array);
                     fortprintf(fd, "         if ((%s %% %s %% ioinfo %% input .and. input_obj %% stream == STREAM_INPUT) .or. &\n", struct_deref, var_ptr->var_array);
                     fortprintf(fd, "             (%s %% %s %% ioinfo %% restart .and. input_obj %% stream == STREAM_RESTART) .or. &\n", struct_deref, var_ptr->var_array);
                     fortprintf(fd, "             (%s %% %s %% ioinfo %% sfc .and. input_obj %% stream == STREAM_SFC)) then\n", struct_deref, var_ptr->var_array);
                     memcpy(var_array, var_ptr->var_array, 1024);
/*                     fortprintf(fd, "            write(0,*) \'exchange halo for %s\'\n", var_ptr->var_array); */
                     fortprintf(fd, "            call mpas_dmpar_exch_halo_field(%s %% %s)\n", struct_deref, var_ptr->var_array);
                     while (var_list_ptr && strncmp(var_array, var_list_ptr->var->var_array, 1024) == 0) {
						var_list_ptr2 = var_list_ptr;
                        var_list_ptr = var_list_ptr->next;
                     }
					 var_list_ptr = var_list_ptr2;
                  }
                  else {
                     fortprintf(fd, "      if (%s %% %s %% isPersistent .and. %s %% %s %% isActive) then\n", struct_deref, var_ptr->name_in_code, struct_deref, var_ptr->name_in_code);
                     fortprintf(fd, "         if ((%s %% %s %% ioinfo %% input .and. input_obj %% stream == STREAM_INPUT) .or. &\n", struct_deref, var_ptr->name_in_code);
                     fortprintf(fd, "             (%s %% %s %% ioinfo %% restart .and. input_obj %% stream == STREAM_RESTART) .or. &\n", struct_deref, var_ptr->name_in_code);
                     fortprintf(fd, "             (%s %% %s %% ioinfo %% sfc .and. input_obj %% stream == STREAM_SFC)) then\n", struct_deref, var_ptr->name_in_code);
/*                     fortprintf(fd, "            write(0,*) \'exchange halo for %s\'\n", var_ptr->name_in_code); */
                     fortprintf(fd, "            call mpas_dmpar_exch_halo_field(%s %% %s)\n", struct_deref, var_ptr->name_in_code);
                  }
            
                  fortprintf(fd, "         end if\n\n");
                  fortprintf(fd, "      end if\n\n");
   
               } else {
                  fortprintf(fd2, "      if (%s %% %s %% isPersistent .and. %s %% %s %% isActive) then\n", struct_deref, var_ptr->name_in_code, struct_deref, var_ptr->name_in_code);
                  fortprintf(fd2, "         if ((%s %% %s %% ioinfo %% input .and. input_obj %% stream == STREAM_INPUT) .or. &\n", struct_deref, var_ptr->name_in_code);
                  fortprintf(fd2, "             (%s %% %s %% ioinfo %% restart .and. input_obj %% stream == STREAM_RESTART) .or. &\n", struct_deref, var_ptr->name_in_code);
                  fortprintf(fd2, "             (%s %% %s %% ioinfo %% sfc .and. input_obj %% stream == STREAM_SFC)) then\n", struct_deref, var_ptr->name_in_code);
				  fortprintf(fd2, "             call mpas_dmpar_copy_field(%s %% %s)\n", struct_deref, var_ptr->name_in_code);
                  fortprintf(fd2, "         end if\n\n");
                  fortprintf(fd2, "      end if\n\n");
			   }
            }
   
            i++;
            dimlist_ptr = dimlist_ptr -> next;
         }
		 }

         if (var_list_ptr) var_list_ptr = var_list_ptr->next;
      }
      group_ptr = group_ptr->next;
   }

   fclose(fd);
   fclose(fd2);

}

void gen_packages(struct package * pkgs){
	FILE * fd;
	struct package * pkg_ptr;

	fd = fopen("define_packages.inc", "w");

	for (pkg_ptr = pkgs; pkg_ptr; pkg_ptr = pkg_ptr->next) {
		if (strlen(pkg_ptr->name) > 0) { 
			fortprintf(fd, "         logical :: %sActive = .false.\n", pkg_ptr->name);
		}
	}

	fclose(fd);
}


void gen_writes(struct group_list * groups, struct variable * vars, struct dimension * dims, struct namelist * namelists)
{
   struct variable * var_ptr;
   struct variable_list * var_list_ptr, *var_list_ptr2;
   struct dimension * dim_ptr;
   struct dimension_list * dimlist_ptr, * lastdim;
   struct group_list * group_ptr;
   struct dtable * dictionary;
   struct namelist * nl;
   FILE * fd;
   char vtype[5];
   char fname[32];
   char struct_deref[1024];
   char var_array[1024];
   char * cp1, * cp2;
   int i, j;
   int ivtype;
   
   
   fd = fopen("add_output_fields.inc", "w");

   group_ptr = groups;
   while (group_ptr) {
      var_list_ptr = group_ptr->vlist;
      while (var_list_ptr) {
         var_ptr = var_list_ptr->var;

         if (group_ptr->ntime_levs > 1)
            snprintf(struct_deref, 1024, "domain %% blocklist %% %s %% time_levs(1) %% %s", group_ptr->name, group_ptr->name);
         else
            snprintf(struct_deref, 1024, "domain %% blocklist %% %s", group_ptr->name);
         
         if (strncmp(var_ptr->var_array, "-", 1024) != 0) {
            fortprintf(fd, "      if ((%s %% %s %% ioinfo %% output .and. output_obj %% stream == OUTPUT) .or. &\n", struct_deref, var_ptr->var_array);
            fortprintf(fd, "          (%s %% %s %% ioinfo %% restart .and. output_obj %% stream == RESTART) .or. &\n", struct_deref, var_ptr->var_array);
            fortprintf(fd, "          (%s %% %s %% ioinfo %% sfc .and. output_obj %% stream == SFC)) then\n", struct_deref, var_ptr->var_array);
            memcpy(var_array, var_ptr->var_array, 1024);
            fortprintf(fd, "         call MPAS_streamAddField(output_obj %% io_stream, %s %% %s, ierr)\n", struct_deref, var_array);
            while (var_list_ptr && strncmp(var_array, var_list_ptr->var->var_array, 1024) == 0) {
			   var_list_ptr2 = var_list_ptr;
               var_list_ptr = var_list_ptr->next;
            }
			var_list_ptr = var_list_ptr2;
         }
         else {
            fortprintf(fd, "      if ((%s %% %s %% ioinfo %% output .and. output_obj %% stream == OUTPUT) .or. &\n", struct_deref, var_ptr->name_in_code);
            fortprintf(fd, "          (%s %% %s %% ioinfo %% restart .and. output_obj %% stream == RESTART) .or. &\n", struct_deref, var_ptr->name_in_code);
            fortprintf(fd, "          (%s %% %s %% ioinfo %% sfc .and. output_obj %% stream == SFC)) then\n", struct_deref, var_ptr->name_in_code);
            fortprintf(fd, "         call MPAS_streamAddField(output_obj %% io_stream, %s %% %s, ierr)\n", struct_deref, var_ptr->name_in_code);
         }
   
         fortprintf(fd, "      end if\n\n");
   
         if (var_list_ptr) var_list_ptr = var_list_ptr->next;
      }
      group_ptr = group_ptr->next;
   }

   fclose(fd);


   fd = fopen("add_output_atts.inc", "w");

   nl = namelists;
   while (nl) {
      if (nl->vtype == LOGICAL) {
         fortprintf(fd, "      if (%s) then\n", nl->name);
         fortprintf(fd, "         call MPAS_writeStreamAtt(output_obj %% io_stream, \'%s\', 'T', ierr)\n", nl->name);
         fortprintf(fd, "      else\n");
         fortprintf(fd, "         call MPAS_writeStreamAtt(output_obj %% io_stream, \'%s\', 'F', ierr)\n", nl->name);
         fortprintf(fd, "      end if\n");
      }
      else {
         fortprintf(fd, "      call MPAS_writeStreamAtt(output_obj %% io_stream, \'%s\', %s, ierr)\n", nl->name, nl->name);
      }
      nl = nl->next;
   }

   fclose(fd);

}
