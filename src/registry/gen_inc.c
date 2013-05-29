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
   char super_array[1024];
   char array_class[1024];
   char outer_dim[1024];
   int i;
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
         memcpy(super_array, var_list_ptr->var->super_array, 1024);
         i = 1;
         while (var_list_ptr) {
            if (strncmp(super_array, var_list_ptr->var->super_array, 1024) != 0) {
               memcpy(super_array, var_list_ptr->var->super_array, 1024);
               i = 1;
             }
            if (strncmp(var_list_ptr->var->array_class, "-", 1024) != 0) fortprintf(fd, "      integer :: index_%s = %i\n", var_list_ptr->var->name_in_code, i++);
            var_list_ptr = var_list_ptr->next;
         }

         var_list_ptr = group_ptr->vlist;
         memcpy(super_array, var_list_ptr->var->super_array, 1024);
         memcpy(array_class, var_list_ptr->var->array_class, 1024);
         class_start = 1;
         class_end = 1;
         i = 1;
         while (var_list_ptr) {
            if (strncmp(var_list_ptr->var->super_array, "-", 1024) != 0) {
               if (strncmp(super_array, var_list_ptr->var->super_array, 1024) != 0) {
                  if (strncmp(super_array, "-", 1024) != 0) fortprintf(fd, "      integer :: %s_end = %i\n", array_class, class_end);
                  if (strncmp(super_array, "-", 1024) != 0) fortprintf(fd, "      integer :: num_%s = %i\n", super_array, i);
                  class_start = 1;
                  class_end = 1;
                  i = 1;
                  memcpy(super_array, var_list_ptr->var->super_array, 1024);
                  memcpy(array_class, var_list_ptr->var->array_class, 1024);
                  fortprintf(fd, "      integer :: %s_start = %i\n", array_class, class_start);
               }
               else if (strncmp(array_class, var_list_ptr->var->array_class, 1024) != 0) {
                  fortprintf(fd, "      integer :: %s_end = %i\n", array_class, class_end);
                  class_start = class_end+1;
                  class_end = class_start;
                  memcpy(array_class, var_list_ptr->var->array_class, 1024);
                  fortprintf(fd, "      integer :: %s_start = %i\n", array_class, class_start);
                  i++;
               }
               else {
                  class_end++;
                  i++;
               }
            }
            var_list_ptr = var_list_ptr->next;
         }
         if (strncmp(super_array, "-", 1024) != 0) fortprintf(fd, "      integer :: %s_end = %i\n", array_class, class_end);
         if (strncmp(super_array, "-", 1024) != 0) fortprintf(fd, "      integer :: num_%s = %i\n", super_array, i);

         var_list_ptr = group_ptr->vlist;
         while (var_list_ptr) {
            var_ptr = var_list_ptr->var;
            if (!strncmp(var_ptr->super_array, "-", 1024)) {
              if (var_ptr->vtype == INTEGER)   fortprintf(fd, "      type (field%idInteger), pointer :: %s\n", var_ptr->ndims, var_ptr->name_in_code);
              if (var_ptr->vtype == REAL)      fortprintf(fd, "      type (field%idReal), pointer :: %s\n", var_ptr->ndims, var_ptr->name_in_code);
              if (var_ptr->vtype == CHARACTER) fortprintf(fd, "      type (field%idChar), pointer :: %s\n", var_ptr->ndims, var_ptr->name_in_code);
            }
            else {
              if (var_ptr->vtype == INTEGER)   fortprintf(fd, "      type (field%idInteger), pointer :: %s\n", var_ptr->ndims+1, var_ptr->super_array);
              if (var_ptr->vtype == REAL)      fortprintf(fd, "      type (field%idReal), pointer :: %s\n", var_ptr->ndims+1, var_ptr->super_array);
              if (var_ptr->vtype == CHARACTER) fortprintf(fd, "      type (field%idChar), pointer :: %s\n", var_ptr->ndims+1, var_ptr->super_array);
              while (var_list_ptr->next && !strncmp(var_list_ptr->next->var->super_array, var_list_ptr->var->super_array, 1024)) var_list_ptr = var_list_ptr->next;
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
         memcpy(super_array, var_list_ptr->var->super_array, 1024);
         i = 1;
         while (var_list_ptr) {
            if (strncmp(super_array, var_list_ptr->var->super_array, 1024) != 0) {
               memcpy(super_array, var_list_ptr->var->super_array, 1024);
               i = 1;
             }
            if (strncmp(var_list_ptr->var->array_class, "-", 1024) != 0) fortprintf(fd, "      integer :: index_%s = %i\n", var_list_ptr->var->name_in_code, i++);
            var_list_ptr = var_list_ptr->next;
         }

         var_list_ptr = group_ptr->vlist;
         sprintf(super_array, "-");
         sprintf(array_class, "-");
         class_start = 1;
         class_end = 1;
         i = 1;

         while (var_list_ptr) {

            /* Is the current variable in a super array? */
            if (strncmp(var_list_ptr->var->super_array, "-", 1024) != 0) {

               /* Have we hit the beginning of a new super array? */
               if (strncmp(super_array, var_list_ptr->var->super_array, 1024) != 0) {
                  /* Finish off the previous super array? */
                  if (strncmp(super_array, "-", 1024) != 0) {
                     fortprintf(fd, "      integer :: %s_end = %i\n", array_class, class_end);
                     fortprintf(fd, "      integer :: num_%s = %i\n", super_array, i);
                  }
                  class_start = 1;
                  class_end = 1;
                  i = 1;
                  memcpy(super_array, var_list_ptr->var->super_array, 1024);
                  memcpy(array_class, var_list_ptr->var->array_class, 1024);
                  fortprintf(fd, "      integer :: %s_start = %i\n", array_class, class_start);
               }
               /* Or have we hit the beginning of a new array class? */
               else if (strncmp(array_class, var_list_ptr->var->array_class, 1024) != 0) {
                  fortprintf(fd, "      integer :: %s_end = %i\n", array_class, class_end);
                  class_start = class_end+1;
                  class_end = class_start;
                  memcpy(array_class, var_list_ptr->var->array_class, 1024);
                  fortprintf(fd, "      integer :: %s_start = %i\n", array_class, class_start);
                  i++;
               }
               else {
                  class_end++;
                  i++;
               }

            }
            var_list_ptr = var_list_ptr->next;

         }
         if (strncmp(super_array, "-", 1024) != 0) fortprintf(fd, "      integer :: %s_end = %i\n", array_class, class_end);
         if (strncmp(super_array, "-", 1024) != 0) fortprintf(fd, "      integer :: num_%s = %i\n", super_array, i);

         var_list_ptr = group_ptr->vlist;
         while (var_list_ptr) {
            var_ptr = var_list_ptr->var;
            if (!strncmp(var_ptr->super_array, "-", 1024)) {
              if (var_ptr->vtype == INTEGER)   fortprintf(fd, "      type (field%idInteger), pointer :: %s\n", var_ptr->ndims, var_ptr->name_in_code);
              if (var_ptr->vtype == REAL)      fortprintf(fd, "      type (field%idReal), pointer :: %s\n", var_ptr->ndims, var_ptr->name_in_code);
              if (var_ptr->vtype == CHARACTER) fortprintf(fd, "      type (field%idChar), pointer :: %s\n", var_ptr->ndims, var_ptr->name_in_code);
            }
            else {
              if (var_ptr->vtype == INTEGER)   fortprintf(fd, "      type (field%idInteger), pointer :: %s\n", var_ptr->ndims+1, var_ptr->super_array);
              if (var_ptr->vtype == REAL)      fortprintf(fd, "      type (field%idReal), pointer :: %s\n", var_ptr->ndims+1, var_ptr->super_array);
              if (var_ptr->vtype == CHARACTER) fortprintf(fd, "      type (field%idChar), pointer :: %s\n", var_ptr->ndims+1, var_ptr->super_array);
              while (var_list_ptr->next && !strncmp(var_list_ptr->next->var->super_array, var_list_ptr->var->super_array, 1024)) var_list_ptr = var_list_ptr->next;
            }
            var_list_ptr = var_list_ptr->next;
         }
   
         fortprintf(fd, "   end type %s_type\n\n\n", group_ptr->name);
   
         if (group_ptr->vlist->var->ntime_levs > 1) {
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
      if (group_ptr->vlist->var->ntime_levs > 1) {
         fortprintf(fd, "      type (%s_multilevel_type), pointer :: %s\n", group_ptr->name, group_ptr->name);
         fortprintf(fd, "      type (%s_type), pointer :: provis\n", group_ptr->name, group_ptr->name);
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
      if (group_ptr->vlist->var->ntime_levs > 1) {
		 fortprintf(fd, "   subroutine mpas_setup_provis_%ss(b)!{{{\n", group_ptr->name);
		 fortprintf(fd, "      type (block_type), pointer :: b\n");
		 fortprintf(fd, "      type (block_type), pointer :: block\n\n");
		 fortprintf(fd, "#include \"dim_dummy_decls_noinput.inc\"\n\n");
		 fortprintf(fd, "      block => b\n");
		 fortprintf(fd, "      do while(associated(block))\n");
		 fortprintf(fd, "#include \"dim_dummy_assigns.inc\"\n\n");
         fortprintf(fd, "         allocate(block %% provis)\n");
         fortprintf(fd, "         call mpas_allocate_%s(block, block %% provis, &\n", group_ptr->name);
         fortprintf(fd, "#include \"dim_dummy_args.inc\"\n");
         fortprintf(fd, "                              )\n\n");
		 fortprintf(fd, "         block => block %% next \n");
		 fortprintf(fd, "      end do\n\n");
		 fortprintf(fd, "      block => b\n");
		 fortprintf(fd, "      do while(associated(block))\n");
         fortprintf(fd, "         if(associated(block %% prev) .and. associated(block %% next)) then\n");
         fortprintf(fd, "            call mpas_create_%s_links(block %% provis, prev = block %% prev %% provis, next = block %% next %% provis)\n", group_ptr->name);
         fortprintf(fd, "         else if(associated(block %% prev)) then\n");
         fortprintf(fd, "            call mpas_create_%s_links(block %% provis, prev = block %% prev %% provis)\n", group_ptr->name);
         fortprintf(fd, "         else if(associated(block %% next)) then\n");
         fortprintf(fd, "            call mpas_create_%s_links(block %% provis, next = block %% next %% provis)\n", group_ptr->name);
         fortprintf(fd, "         else\n");
         fortprintf(fd, "            call mpas_create_%s_links(block %% provis)\n", group_ptr->name);
         fortprintf(fd, "         end if\n");
		 fortprintf(fd, "         block => block %% next \n");
		 fortprintf(fd, "      end do\n");
		 fortprintf(fd, "   end subroutine mpas_setup_provis_%ss!}}}\n\n", group_ptr->name);

		 fortprintf(fd, "   subroutine mpas_deallocate_provis_%ss(b)!{{{\n", group_ptr->name);
		 fortprintf(fd, "      type (block_type), pointer :: b\n");
		 fortprintf(fd, "      type (block_type), pointer :: block\n\n");
		 fortprintf(fd, "      block => b\n");
		 fortprintf(fd, "      do while(associated(block))\n");
		 fortprintf(fd, "         call mpas_deallocate_%s(block %% provis)\n", group_ptr->name);
		 fortprintf(fd, "         deallocate(block %% provis)\n");
		 fortprintf(fd, "         block => block %% next\n");
		 fortprintf(fd, "      end do\n");
		 fortprintf(fd, "   end subroutine mpas_deallocate_provis_%ss!}}}\n", group_ptr->name);
	  }
      group_ptr = group_ptr->next;
   }
   fclose(fd);



   /* To be included in allocate_block */
   fd = fopen("block_allocs.inc", "w");
   group_ptr = groups;
   while (group_ptr) {
      fortprintf(fd, "      allocate(b %% %s)\n", group_ptr->name);
      if (group_ptr->vlist->var->ntime_levs > 1) {
         fortprintf(fd, "      b %% %s %% nTimeLevels = %i\n", group_ptr->name, group_ptr->vlist->var->ntime_levs);
         fortprintf(fd, "      allocate(b %% %s %% time_levs(%i))\n", group_ptr->name, group_ptr->vlist->var->ntime_levs);
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
      if (group_ptr->vlist->var->ntime_levs > 1) {
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
         if (strncmp(var_ptr->super_array, "-", 1024) != 0) {
            memcpy(super_array, var_ptr->super_array, 1024);
            memcpy(array_class, var_ptr->array_class, 1024);
            i = 0;
            while (var_list_ptr && strncmp(super_array, var_list_ptr->var->super_array, 1024) == 0) {
               i++;
               var_list_ptr2 = var_list_ptr;
               var_list_ptr = var_list_ptr->next;
            }
            var_ptr2 = var_list_ptr2->var;
            fortprintf(fd, "      allocate(%s %% %s)\n", group_ptr->name, var_ptr2->super_array);
            fortprintf(fd, "      allocate(%s %% %s %% ioinfo)\n", group_ptr->name, var_ptr2->super_array);
            fortprintf(fd, "      %s %% %s %% fieldName = \'%s\'\n", group_ptr->name, var_ptr2->super_array, var_ptr2->super_array);
            fortprintf(fd, "      %s %% %s %% isSuperArray = .true.\n", group_ptr->name, var_ptr2->super_array);
            fortprintf(fd, "      allocate(%s %% %s %% constituentNames(%i))\n", group_ptr->name, var_ptr2->super_array, i);

            /* Initialization for constituent names */
            i = 0;
            var_list_ptr3 = group_ptr->vlist;
            while (var_list_ptr3) {
               if (strncmp(super_array, var_list_ptr3->var->super_array, 1024) == 0) {
                  i++;
                  fortprintf(fd, "      %s %% %s %% constituentNames(%i) = \'%s\'\n", group_ptr->name, var_ptr2->super_array, i, var_list_ptr3->var->name_in_file);
               }
               var_list_ptr3 = var_list_ptr3->next;
            }

			if(var_ptr2->persistence == PERSISTENT){
               fortprintf(fd, "      allocate(%s %% %s %% array(%i, ", group_ptr->name, var_ptr2->super_array, i);
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
               if (var_ptr->vtype == INTEGER)
                  fortprintf(fd, "      %s %% %s %% array = 0\n", group_ptr->name, var_ptr2->super_array ); /* initialize field to zero */
               else if (var_ptr->vtype == REAL)
                  fortprintf(fd, "      %s %% %s %% array = 0.0\n", group_ptr->name, var_ptr2->super_array ); /* initialize field to zero */
               else if (var_ptr->vtype == CHARACTER)
                  fortprintf(fd, "      %s %% %s %% array = \'\'\n", group_ptr->name, var_ptr2->super_array ); /* initialize field to zero */
			}

            fortprintf(fd, "      %s %% %s %% dimSizes(1) = %i\n", group_ptr->name, var_ptr2->super_array, i);
            fortprintf(fd, "      %s %% %s %% dimNames(1) = \'num_%s\'\n", group_ptr->name, var_ptr2->super_array, var_ptr2->super_array);
            dimlist_ptr = var_ptr2->dimlist;
            i = 2;
            while (dimlist_ptr) {
               if (!strncmp(dimlist_ptr->dim->name_in_file, "nCells", 1024) ||
                   !strncmp(dimlist_ptr->dim->name_in_file, "nEdges", 1024) ||
                   !strncmp(dimlist_ptr->dim->name_in_file, "nVertices", 1024))
                  if (!dimlist_ptr->dim->namelist_defined) {
					 if (var_ptr2->persistence == PERSISTENT){
                        fortprintf(fd, "      %s %% %s %% dimSizes(%i) = %s\n", group_ptr->name, var_ptr2->super_array, i, dimlist_ptr->dim->name_in_code);
                        fortprintf(fd, "      %s %% %s %% dimNames(%i) = \'%s\'\n", group_ptr->name, var_ptr2->super_array, i, dimlist_ptr->dim->name_in_file);
					 } 
					 else {
                        fortprintf(fd, "      %s %% %s %% dimSizes(%i) = %s+1\n", group_ptr->name, var_ptr2->super_array, i, dimlist_ptr->dim->name_in_code);
                        fortprintf(fd, "      %s %% %s %% dimNames(%i) = \'%s\'\n", group_ptr->name, var_ptr2->super_array, i, dimlist_ptr->dim->name_in_file);
					 }
                  }
                  else {
                     fortprintf(fd, "      %s %% %s %% dimSizes(%i) = %s\n", group_ptr->name, var_ptr2->super_array, i, dimlist_ptr->dim->name_in_file);
                     fortprintf(fd, "      %s %% %s %% dimNames(%i) = \'%s\'\n", group_ptr->name, var_ptr2->super_array, i, dimlist_ptr->dim->name_in_file);
                  }
               else
                  if (dimlist_ptr->dim->namelist_defined) {
                     fortprintf(fd, "      %s %% %s %% dimSizes(%i) = %s\n", group_ptr->name, var_ptr2->super_array, i, dimlist_ptr->dim->name_in_file);
                     fortprintf(fd, "      %s %% %s %% dimNames(%i) = \'%s\'\n", group_ptr->name, var_ptr2->super_array, i, dimlist_ptr->dim->name_in_file);
                  }
                  else {
                     fortprintf(fd, "      %s %% %s %% dimSizes(%i) = %s\n", group_ptr->name, var_ptr2->super_array, i, dimlist_ptr->dim->name_in_code);
                     fortprintf(fd, "      %s %% %s %% dimNames(%i) = \'%s\'\n", group_ptr->name, var_ptr2->super_array, i, dimlist_ptr->dim->name_in_file);
                  }
               i++;
               dimlist_ptr = dimlist_ptr->next;
            }
            if (var_ptr2->timedim) fortprintf(fd, "      %s %% %s %% hasTimeDimension = .true.\n", group_ptr->name, var_ptr2->super_array);
            else fortprintf(fd, "      %s %% %s %% hasTimeDimension = .false.\n", group_ptr->name, var_ptr2->super_array);
            fortprintf(fd, "      nullify(%s %% %s %% prev)\n", group_ptr->name, var_ptr2->super_array);
            fortprintf(fd, "      nullify(%s %% %s %% next)\n", group_ptr->name, var_ptr2->super_array);
            fortprintf(fd, "      nullify(%s %% %s %% sendList)\n", group_ptr->name, var_ptr2->super_array);
            fortprintf(fd, "      nullify(%s %% %s %% recvList)\n", group_ptr->name, var_ptr2->super_array);
            fortprintf(fd, "      nullify(%s %% %s %% copyList)\n", group_ptr->name, var_ptr2->super_array);

            if (var_ptr2->iostreams & INPUT0) 
               fortprintf(fd, "      %s %% %s %% ioinfo %% input = .true.\n", group_ptr->name, var_ptr2->super_array);
            else
               fortprintf(fd, "      %s %% %s %% ioinfo %% input = .false.\n", group_ptr->name, var_ptr2->super_array);

            if (var_ptr2->iostreams & SFC0) 
               fortprintf(fd, "      %s %% %s %% ioinfo %% sfc = .true.\n", group_ptr->name, var_ptr2->super_array);
            else
               fortprintf(fd, "      %s %% %s %% ioinfo %% sfc = .false.\n", group_ptr->name, var_ptr2->super_array);

            if (var_ptr2->iostreams & RESTART0) 
               fortprintf(fd, "      %s %% %s %% ioinfo %% restart = .true.\n", group_ptr->name, var_ptr2->super_array);
            else
               fortprintf(fd, "      %s %% %s %% ioinfo %% restart = .false.\n", group_ptr->name, var_ptr2->super_array);

            if (var_ptr2->iostreams & OUTPUT0) 
               fortprintf(fd, "      %s %% %s %% ioinfo %% output = .true.\n", group_ptr->name, var_ptr2->super_array);
            else
               fortprintf(fd, "      %s %% %s %% ioinfo %% output = .false.\n", group_ptr->name, var_ptr2->super_array);

            fortprintf(fd, "      %s %% %s %% block => b\n", group_ptr->name, var_ptr2->super_array);
            fortprintf(fd, "\n");
         }
         else {
            fortprintf(fd, "      allocate(%s %% %s)\n", group_ptr->name, var_ptr->name_in_code);
            fortprintf(fd, "      allocate(%s %% %s %% ioinfo)\n", group_ptr->name, var_ptr->name_in_code);
            fortprintf(fd, "      %s %% %s %% fieldName = \'%s\'\n", group_ptr->name, var_ptr->name_in_code, var_ptr->name_in_file);
            fortprintf(fd, "      %s %% %s %% isSuperArray = .false.\n", group_ptr->name, var_ptr->name_in_code);
            if (var_ptr->ndims > 0) {
	  		  if(var_ptr->persistence == SCRATCH){
				  fortprintf(fd, "      %s %% %s %% isPersistent = .false.\n", group_ptr->name, var_ptr->name_in_code);
				  fortprintf(fd, "      nullify(%s %% %s %% array)\n", group_ptr->name, var_ptr->name_in_code); 
			  } else if(var_ptr->persistence == PERSISTENT){
				  fortprintf(fd, "      %s %% %s %% isPersistent = .true.\n", group_ptr->name, var_ptr->name_in_code);
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
               if (var_ptr->vtype == INTEGER)
                  fortprintf(fd, "      %s %% %s %% array = 0\n", group_ptr->name, var_ptr->name_in_code ); /* initialize field to zero */
               else if (var_ptr->vtype == REAL)
                  fortprintf(fd, "      %s %% %s %% array = 0.0\n", group_ptr->name, var_ptr->name_in_code ); /* initialize field to zero */
               else if (var_ptr->vtype == CHARACTER)
                  fortprintf(fd, "      %s %% %s %% array = \'\'\n", group_ptr->name, var_ptr->name_in_code ); /* initialize field to zero */

			  }
               dimlist_ptr = var_ptr->dimlist;
               i = 1;
               while (dimlist_ptr) {
                  if (!strncmp(dimlist_ptr->dim->name_in_file, "nCells", 1024) ||
                      !strncmp(dimlist_ptr->dim->name_in_file, "nEdges", 1024) ||
                      !strncmp(dimlist_ptr->dim->name_in_file, "nVertices", 1024))
                     if (!dimlist_ptr->dim->namelist_defined) {
						if(var_ptr->persistence == PERSISTENT){
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
         if (strncmp(var_ptr->super_array, "-", 1024) != 0) {
            memcpy(super_array, var_ptr->super_array, 1024);
            memcpy(array_class, var_ptr->array_class, 1024);
            i = 0;
            while (var_list_ptr && strncmp(super_array, var_list_ptr->var->super_array, 1024) == 0) {
               i++;
               var_list_ptr2 = var_list_ptr;
               var_list_ptr = var_list_ptr->next;
            }
            fortprintf(fd, "      if(associated(%s %% %s %% array)) then\n", group_ptr->name, var_list_ptr2->var->super_array);
            fortprintf(fd, "         deallocate(%s %% %s %% array)\n", group_ptr->name, var_list_ptr2->var->super_array);
            fortprintf(fd, "      end if\n");
            fortprintf(fd, "      deallocate(%s %% %s %% ioinfo)\n", group_ptr->name, var_list_ptr2->var->super_array);
            fortprintf(fd, "      call mpas_deallocate_attlist(%s %% %s %% attList)\n", group_ptr->name, var_list_ptr2->var->super_array);
            fortprintf(fd, "      deallocate(%s %% %s)\n\n", group_ptr->name, var_list_ptr2->var->super_array);
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
         if (strncmp(var_ptr->super_array, "-", 1024) != 0) {
            memcpy(super_array, var_ptr->super_array, 1024);
            memcpy(array_class, var_ptr->array_class, 1024);
            i = 0;
            while (var_list_ptr && strncmp(super_array, var_list_ptr->var->super_array, 1024) == 0) {
               i++;
               var_list_ptr2 = var_list_ptr;
               var_list_ptr = var_list_ptr->next;
            }
            var_ptr2 = var_list_ptr2->var;
            if (var_ptr2->ndims > 0) 
               fortprintf(fd, "      dest %% %s %% array = src %% %s %% array\n", var_ptr2->super_array, var_ptr2->super_array);
            else
               fortprintf(fd, "      dest %% %s %% scalar = src %% %s %% scalar\n", var_ptr2->super_array, var_ptr2->super_array);
         }
         else {
			if (var_ptr->persistence == PERSISTENT){
            if (var_ptr->ndims > 0) 
               fortprintf(fd, "      dest %% %s %% array = src %% %s %% array\n", var_ptr->name_in_code, var_ptr->name_in_code);
            else
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
      if (group_ptr->vlist->var->ntime_levs > 1) {
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

            if (strncmp(var_ptr->super_array, "-", 1024) != 0) 
            {
               if (var_ptr->vtype == INTEGER) sprintf(type_str, "int%id", var_ptr->ndims+1); 
               else if (var_ptr->vtype == REAL) sprintf(type_str, "real%id", var_ptr->ndims+1); 
               else if (var_ptr->vtype == CHARACTER) sprintf(type_str, "char%id", var_ptr->ndims+1); 

               memcpy(super_array, var_ptr->super_array, 1024);

               while (var_list_ptr && strncmp(super_array, var_list_ptr->var->super_array, 1024) == 0)
               {
                  var_list_ptr2 = var_list_ptr;
                  var_list_ptr = var_list_ptr->next;
               }
               var_ptr2 = var_list_ptr2->var;

               fortprintf(fd, "      %s => %s %% time_levs(1) %% %s %% %s %% array\n", type_str, group_ptr->name, group_ptr->name, var_ptr2->super_array);

               fortprintf(fd, "      do i=1,%s %% nTimeLevels-1\n", group_ptr->name);
               fortprintf(fd, "         %s %% time_levs(i) %% %s %% %s %% array => %s %% time_levs(i+1) %% %s %% %s %% array\n", group_ptr->name, group_ptr->name, var_ptr2->super_array, group_ptr->name, group_ptr->name, var_ptr2->super_array);
               fortprintf(fd, "      end do\n");

               fortprintf(fd, "      %s %% time_levs(%s %% nTimeLevels) %% %s %% %s %% array=> %s\n\n", group_ptr->name, group_ptr->name, group_ptr->name, var_ptr2->super_array, type_str);
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
   fortprintf(fd, "      subroutine mpas_create_field_links(b)\n\n");
   fortprintf(fd, "         implicit none\n");
   fortprintf(fd, "         type (block_type), pointer :: b\n");
   fortprintf(fd, "         type (block_type), pointer :: prev, next\n\n");
   fortprintf(fd, "         if(associated(b %% prev)) then\n");
   fortprintf(fd, "           prev => b %% prev\n");
   fortprintf(fd, "         else\n");
   fortprintf(fd, "           nullify(prev)\n");
   fortprintf(fd, "         end if\n");
   fortprintf(fd, "         if(associated(b %% next)) then\n");
   fortprintf(fd, "           next => b %% next\n");
   fortprintf(fd, "         else\n");
   fortprintf(fd, "           nullify(next)\n");
   fortprintf(fd, "         end if\n\n");
   group_ptr = groups;
   while (group_ptr)
   {
     var_list_ptr = group_ptr->vlist;
     var_list_ptr = var_list_ptr->next;

     if (!var_list_ptr) break;

     var_ptr = var_list_ptr->var;
     
     int ntime_levs = 1;
     
     if (strncmp(var_ptr->super_array, "-", 1024) != 0) 
     {
         memcpy(super_array, var_ptr->super_array, 1024);
         memcpy(array_class, var_ptr->array_class, 1024);
         while (var_list_ptr && strncmp(super_array, var_list_ptr->var->super_array, 1024) == 0)
         {
            var_list_ptr2 = var_list_ptr;
            var_list_ptr = var_list_ptr->next;
         }
         var_ptr2 = var_list_ptr2->var;
         get_outer_dim(var_ptr2, outer_dim);
         ntime_levs = var_ptr2->ntime_levs;

         if(ntime_levs > 1)
         {
            for(i=1; i<=ntime_levs; i++) 
            {
				fortprintf(fd, "         if(associated(next) .and. associated(prev)) then\n");	
//				fortprintf(fd, "           call mpas_create_%s_links(b %% %s %% time_levs(%i) %% %s, prev = prev %% %s %% time_levs(%i) %% %s, next = next %% %s %% time_levs(%i) %% %s)\n", group_ptr->name, group_ptr->name, i, group_ptr->name, i, group_ptr->name, group_ptr->name, i, group_ptr->name);
				fortprintf(fd, "           call mpas_create_%s_links(b %% %s %% time_levs(%i) %% %s, ", group_ptr->name, group_ptr->name, i, group_ptr->name, i);
				fortprintf(fd, " prev = prev %% %s %% time_levs(%i) %% %s,", group_ptr->name, i, group_ptr->name);
				fortprintf(fd, " next = next %% %s %% time_levs(%i) %% %s)\n", group_ptr->name, i, group_ptr->name);
				fortprintf(fd, "         else if(associated(next)) then\n");	
				fortprintf(fd, "           call mpas_create_%s_links(b %% %s %% time_levs(%i) %% %s, next = next %% %s %% time_levs(%i) %% %s)\n", group_ptr->name, group_ptr->name, i, group_ptr->name, group_ptr->name, i, group_ptr->name);
				fortprintf(fd, "         else if(associated(prev)) then\n");	
				fortprintf(fd, "           call mpas_create_%s_links(b %% %s %% time_levs(%i) %% %s, prev = prev %% %s %% time_levs(%i) %% %s)\n", group_ptr->name, group_ptr->name, i, group_ptr->name, group_ptr->name, i, group_ptr->name);
				fortprintf(fd, "         else\n");
				fortprintf(fd, "           call mpas_create_%s_links(b %% %s %% time_levs(%i) %% %s)\n", group_ptr->name, group_ptr->name, i, group_ptr->name);
				fortprintf(fd, "         end if\n\n");
            }	
         }
         else
         {
			fortprintf(fd, "         if(associated(next) .and. associated(prev)) then\n");	
            fortprintf(fd, "           call mpas_create_%s_links(b %% %s, prev = prev %% %s, next = next %% %s)\n", group_ptr->name, group_ptr->name, group_ptr->name, group_ptr->name); 
			fortprintf(fd, "         else if(associated(next)) then\n");	
            fortprintf(fd, "           call mpas_create_%s_links(b %% %s, next = next %% %s)\n", group_ptr->name, group_ptr->name, group_ptr->name); 
			fortprintf(fd, "         else if(associated(prev)) then\n");	
            fortprintf(fd, "           call mpas_create_%s_links(b %% %s, prev = prev %% %s)\n", group_ptr->name, group_ptr->name, group_ptr->name); 
			fortprintf(fd, "         else\n");
            fortprintf(fd, "           call mpas_create_%s_links(b %% %s)\n", group_ptr->name, group_ptr->name); 
			fortprintf(fd, "         end if\n\n");
         }
     }
     else if (var_ptr->ndims > 0)
     {
         get_outer_dim(var_ptr, outer_dim);
         ntime_levs = var_ptr->ntime_levs;

         if(ntime_levs > 1)
         {
            for(i=1; i<=ntime_levs; i++) 
            {
				fortprintf(fd, "         if(associated(next) .and. associated(prev)) then\n");	
				fortprintf(fd, "           call mpas_create_%s_links(b %% %s %% time_levs(%i) %% %s, prev = prev %% %s %% time_levs(%i) %% %s, next = next %% %s %% time_levs(%i) %% %s)\n", group_ptr->name, group_ptr->name, i, group_ptr->name, group_ptr->name, i, group_ptr->name, group_ptr->name, i, group_ptr->name);
				fortprintf(fd, "         else if(associated(next)) then\n");	
				fortprintf(fd, "           call mpas_create_%s_links(b %% %s %% time_levs(%i) %% %s, next = next %% %s %% time_levs(%i) %% %s)\n", group_ptr->name, group_ptr->name, i, group_ptr->name, group_ptr->name, i, group_ptr->name);
				fortprintf(fd, "         else if(associated(prev)) then\n");	
				fortprintf(fd, "           call mpas_create_%s_links(b %% %s %% time_levs(%i) %% %s, prev = prev %% %s %% time_levs(%i) %% %s)\n", group_ptr->name, group_ptr->name, i, group_ptr->name, group_ptr->name, i, group_ptr->name);
				fortprintf(fd, "         else\n");
				fortprintf(fd, "           call mpas_create_%s_links(b %% %s %% time_levs(%i) %% %s)\n", group_ptr->name, group_ptr->name, i, group_ptr->name);
				fortprintf(fd, "         end if\n\n");
            }	
         }
         else
         {
			 fortprintf(fd, "         if(associated(next) .and. associated(prev)) then\n");	
			 fortprintf(fd, "           call mpas_create_%s_links(b %% %s, prev = prev %% %s, next = next %% %s)\n", group_ptr->name, group_ptr->name, group_ptr->name, group_ptr->name);
			 fortprintf(fd, "         else if(associated(next)) then\n");	
			 fortprintf(fd, "           call mpas_create_%s_links(b %% %s, next = next %% %s)\n", group_ptr->name, group_ptr->name, group_ptr->name); 
			 fortprintf(fd, "         else if(associated(prev)) then\n");	
			 fortprintf(fd, "           call mpas_create_%s_links(b %% %s, prev = prev %% %s)\n", group_ptr->name, group_ptr->name, group_ptr->name); 
			 fortprintf(fd, "         else\n");
			 fortprintf(fd, "           call mpas_create_%s_links(b %% %s)\n", group_ptr->name, group_ptr->name); 
			 fortprintf(fd, "         end if\n\n");
         }
     }

     group_ptr = group_ptr->next;
   }
   fortprintf(fd, "\n      end subroutine mpas_create_field_links\n\n\n");

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
         if (strncmp(var_ptr->super_array, "-", 1024) != 0) {
            memcpy(super_array, var_ptr->super_array, 1024);
            memcpy(array_class, var_ptr->array_class, 1024);
            while (var_list_ptr && strncmp(super_array, var_list_ptr->var->super_array, 1024) == 0) {
               var_list_ptr2 = var_list_ptr;
               var_list_ptr = var_list_ptr->next;
            }
            var_ptr2 = var_list_ptr2->var;
            get_outer_dim(var_ptr2, outer_dim);
            
               if (strncmp("nCells",outer_dim,1024) == 0) {
                  fortprintf(fd, "         %s %% %s %% sendList => %s %% %s %% block %% parinfo %% cellsToSend\n", group_ptr->name, var_ptr2->super_array, group_ptr->name, var_ptr2->super_array);
                  fortprintf(fd, "         %s %% %s %% recvList => %s %% %s %% block %% parinfo %% cellsToRecv\n", group_ptr->name, var_ptr2->super_array, group_ptr->name, var_ptr2->super_array);
                  fortprintf(fd, "         %s %% %s %% copyList => %s %% %s %% block %% parinfo %% cellsToCopy\n", group_ptr->name, var_ptr2->super_array, group_ptr->name, var_ptr2->super_array);
				  fortprintf(fd, "         if(present(prev)) then\n");
				  fortprintf(fd, "           %s %% %s %% prev => prev %% %s\n", group_ptr->name, var_ptr2->super_array, var_ptr2->super_array);
				  fortprintf(fd, "         else\n");
				  fortprintf(fd, "           nullify(%s %% %s %% prev)\n", group_ptr->name, var_ptr2->super_array);
				  fortprintf(fd, "         end if\n");
				  fortprintf(fd, "         if(present(next)) then\n");
				  fortprintf(fd, "           %s %% %s %% next => next %% %s\n", group_ptr->name, var_ptr2->super_array, var_ptr2->super_array);
				  fortprintf(fd, "         else\n");
				  fortprintf(fd, "           nullify(%s %% %s %% next)\n", group_ptr->name, var_ptr2->super_array);
				  fortprintf(fd, "         end if\n\n");
               }
               else if (strncmp("nEdges",outer_dim,1024) == 0) {
                  fortprintf(fd, "         %s %% %s %% sendList => %s %% %s %% block %% parinfo %% edgesToSend\n", group_ptr->name, var_ptr2->super_array, group_ptr->name, var_ptr2->super_array);
                  fortprintf(fd, "         %s %% %s %% recvList => %s %% %s %% block %% parinfo %% edgesToRecv\n", group_ptr->name, var_ptr2->super_array, group_ptr->name, var_ptr2->super_array);
                  fortprintf(fd, "         %s %% %s %% copyList => %s %% %s %% block %% parinfo %% edgesToCopy\n", group_ptr->name, var_ptr2->super_array, group_ptr->name, var_ptr2->super_array);
				  fortprintf(fd, "         if(present(prev)) then\n");
				  fortprintf(fd, "           %s %% %s %% prev => prev %% %s\n", group_ptr->name, var_ptr2->super_array, var_ptr2->super_array);
				  fortprintf(fd, "         else\n");
				  fortprintf(fd, "           nullify(%s %% %s %% prev)\n", group_ptr->name, var_ptr2->super_array);
				  fortprintf(fd, "         end if\n");
				  fortprintf(fd, "         if(present(next)) then\n");
				  fortprintf(fd, "           %s %% %s %% next => next %% %s\n", group_ptr->name, var_ptr2->super_array, var_ptr2->super_array);
				  fortprintf(fd, "         else\n");
				  fortprintf(fd, "           nullify(%s %% %s %% next)\n", group_ptr->name, var_ptr2->super_array);
				  fortprintf(fd, "         end if\n\n");
               }
               else if (strncmp("nVertices",outer_dim,1024) == 0) {
                  fortprintf(fd, "         %s %% %s %% sendList => %s %% %s %% block %% parinfo %% verticesToSend\n", group_ptr->name, var_ptr2->super_array, group_ptr->name, var_ptr2->super_array);
                  fortprintf(fd, "         %s %% %s %% recvList => %s %% %s %% block %% parinfo %% verticesToRecv\n", group_ptr->name, var_ptr2->super_array, group_ptr->name, var_ptr2->super_array);
                  fortprintf(fd, "         %s %% %s %% copyList => %s %% %s %% block %% parinfo %% verticesToCopy\n", group_ptr->name, var_ptr2->super_array, group_ptr->name, var_ptr2->super_array);
				  fortprintf(fd, "         if(present(prev)) then\n");
				  fortprintf(fd, "           %s %% %s %% prev => prev %% %s\n", group_ptr->name, var_ptr2->super_array, var_ptr2->super_array);
				  fortprintf(fd, "         else\n");
				  fortprintf(fd, "           nullify(%s %% %s %% prev)\n", group_ptr->name, var_ptr2->super_array);
				  fortprintf(fd, "         end if\n");
				  fortprintf(fd, "         if(present(next)) then\n");
				  fortprintf(fd, "           %s %% %s %% next => next %% %s\n", group_ptr->name, var_ptr2->super_array, var_ptr2->super_array);
				  fortprintf(fd, "         else\n");
				  fortprintf(fd, "           nullify(%s %% %s %% next)\n", group_ptr->name, var_ptr2->super_array);
				  fortprintf(fd, "         end if\n\n");
               } else {
				  fortprintf(fd, "         nullify(%s %% %s %% sendList)\n", group_ptr->name, var_ptr2->super_array);
				  fortprintf(fd, "         nullify(%s %% %s %% recvList)\n", group_ptr->name, var_ptr2->super_array);
				  fortprintf(fd, "         nullify(%s %% %s %% copyList)\n", group_ptr->name, var_ptr2->super_array);
				  fortprintf(fd, "         if(present(prev)) then\n");
				  fortprintf(fd, "           %s %% %s %% prev => prev %% %s\n", group_ptr->name, var_ptr2->super_array, var_ptr2->super_array);
				  fortprintf(fd, "         else\n");
				  fortprintf(fd, "           nullify(%s %% %s %% prev)\n", group_ptr->name, var_ptr2->super_array);
				  fortprintf(fd, "         end if\n");
				  fortprintf(fd, "         if(present(next)) then\n");
				  fortprintf(fd, "           %s %% %s %% next => next %% %s\n", group_ptr->name, var_ptr2->super_array, var_ptr2->super_array);
				  fortprintf(fd, "         else\n");
				  fortprintf(fd, "           nullify(%s %% %s %% next)\n", group_ptr->name, var_ptr2->super_array);
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
   char super_array[1024];
   char struct_deref[1024];
   char * cp1, * cp2;
   int i, j;
   int ivtype;


#ifdef LEGACY_CODE
   /*
    *  Generate declarations of IDs belonging in io_input_object
    */
   fd = fopen("io_input_obj_decls.inc", "w");

   fortprintf(fd, "      integer :: rdDimIDTime\n");
   dim_ptr = dims;
   while (dim_ptr) {
      if (dim_ptr->constant_value < 0 && !dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) fortprintf(fd, "      integer :: rdDimID%s\n", dim_ptr->name_in_file);
      dim_ptr = dim_ptr->next;
   }
   fortprintf(fd, "\n");

   fortprintf(fd, "      integer :: rdLocalTime\n");
   dim_ptr = dims;
   while (dim_ptr) {
      if (dim_ptr->constant_value < 0 && !dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) fortprintf(fd, "      integer :: rdLocal%s\n", dim_ptr->name_in_file);
      dim_ptr = dim_ptr->next;
   }
   fortprintf(fd, "\n");

   var_ptr = vars;
   while (var_ptr) {
      fortprintf(fd, "      integer :: rdVarID%s\n", var_ptr->name_in_file);
      var_ptr = var_ptr->next;
   }

   fclose(fd);


   /*
    *  Definitions of read bounds and exchange lists for non-decomposed fields
    */
   fd = fopen("nondecomp_dims.inc", "w");

   dim_ptr = dims;
   while (dim_ptr) {

      if (strncmp(dim_ptr->name_in_file,"nCells",1024) != 0 &&
          strncmp(dim_ptr->name_in_file,"nEdges",1024) != 0 &&
          strncmp(dim_ptr->name_in_file,"nVertices",1024) != 0 &&
          strncmp(dim_ptr->name_in_file,"nVertLevels",1024) != 0
         ) {

         if (is_derived_dim(dim_ptr->name_in_code)) {
            fortprintf(fd, "      integer :: read%sStart\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "      integer :: read%sCount\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "      type (exchange_list), pointer :: send%sList\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "      type (exchange_list), pointer :: recv%sList\n", dim_ptr->name_in_file+1);
         }
         else if (dim_ptr->constant_value > 0) {
            fortprintf(fd, "      integer :: read%sStart\n", dim_ptr->name_in_file);
            fortprintf(fd, "      integer :: read%sCount\n", dim_ptr->name_in_file);
            fortprintf(fd, "      type (exchange_list), pointer :: send%sList\n", dim_ptr->name_in_file);
            fortprintf(fd, "      type (exchange_list), pointer :: recv%sList\n", dim_ptr->name_in_file);
         }
         else if (dim_ptr->namelist_defined) {
            fortprintf(fd, "      integer :: read%sStart\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "      integer :: read%sCount\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "      type (exchange_list), pointer :: send%sList\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "      type (exchange_list), pointer :: recv%sList\n", dim_ptr->name_in_file+1);
         }
         else {
            fortprintf(fd, "      integer :: read%sStart\n", dim_ptr->name_in_code+1);
            fortprintf(fd, "      integer :: read%sCount\n", dim_ptr->name_in_code+1);
            fortprintf(fd, "      type (exchange_list), pointer :: send%sList\n", dim_ptr->name_in_code+1);
            fortprintf(fd, "      type (exchange_list), pointer :: recv%sList\n", dim_ptr->name_in_code+1);
         }
      }

      dim_ptr = dim_ptr->next;
   }

   fortprintf(fd, "\n");

   dim_ptr = dims;
   while (dim_ptr) {

      if (strncmp(dim_ptr->name_in_file,"nCells",1024) != 0 &&
          strncmp(dim_ptr->name_in_file,"nEdges",1024) != 0 &&
          strncmp(dim_ptr->name_in_file,"nVertices",1024) != 0 &&
          strncmp(dim_ptr->name_in_file,"nVertLevels",1024) != 0
         ) {

         if (is_derived_dim(dim_ptr->name_in_code)) {
            fortprintf(fd, "      read%sStart = 1\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "      read%sCount = block %% mesh %% %s\n", dim_ptr->name_in_file+1, dim_ptr->name_in_code);
            fortprintf(fd, "      allocate(send%sList)\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "      allocate(recv%sList)\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "      nullify(send%sList %% next)\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "      nullify(recv%sList %% next)\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "      recv%sList %% procID = dminfo %% my_proc_id\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "      send%sList %% procID = dminfo %% my_proc_id\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "      recv%sList %% nlist = read%sCount\n", dim_ptr->name_in_file+1, dim_ptr->name_in_file+1);
            fortprintf(fd, "      send%sList %% nlist = read%sCount\n", dim_ptr->name_in_file+1, dim_ptr->name_in_file+1);
            fortprintf(fd, "      allocate(recv%sList %% list(read%sCount))\n", dim_ptr->name_in_file+1, dim_ptr->name_in_file+1);
            fortprintf(fd, "      allocate(send%sList %% list(read%sCount))\n", dim_ptr->name_in_file+1, dim_ptr->name_in_file+1);
            fortprintf(fd, "      do i=1,read%sCount\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "         recv%sList %% list(i) = i\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "         send%sList %% list(i) = i\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "      end do\n");
            fortprintf(fd, "\n");
         }
         else if (dim_ptr->constant_value > 0) {
            fortprintf(fd, "      read%sStart = 1\n", dim_ptr->name_in_file);
            fortprintf(fd, "      read%sCount = %s\n", dim_ptr->name_in_file, dim_ptr->name_in_code);
            fortprintf(fd, "      allocate(send%sList)\n", dim_ptr->name_in_file);
            fortprintf(fd, "      allocate(recv%sList)\n", dim_ptr->name_in_file);
            fortprintf(fd, "      nullify(send%sList %% next)\n", dim_ptr->name_in_file);
            fortprintf(fd, "      nullify(recv%sList %% next)\n", dim_ptr->name_in_file);
            fortprintf(fd, "      recv%sList %% procID = dminfo %% my_proc_id\n", dim_ptr->name_in_file);
            fortprintf(fd, "      send%sList %% procID = dminfo %% my_proc_id\n", dim_ptr->name_in_file);
            fortprintf(fd, "      recv%sList %% nlist = read%sCount\n", dim_ptr->name_in_file, dim_ptr->name_in_file);
            fortprintf(fd, "      send%sList %% nlist = read%sCount\n", dim_ptr->name_in_file, dim_ptr->name_in_file);
            fortprintf(fd, "      allocate(recv%sList %% list(read%sCount))\n", dim_ptr->name_in_file, dim_ptr->name_in_file);
            fortprintf(fd, "      allocate(send%sList %% list(read%sCount))\n", dim_ptr->name_in_file, dim_ptr->name_in_file);
            fortprintf(fd, "      do i=1,read%sCount\n", dim_ptr->name_in_file);
            fortprintf(fd, "         recv%sList %% list(i) = i\n", dim_ptr->name_in_file);
            fortprintf(fd, "         send%sList %% list(i) = i\n", dim_ptr->name_in_file);
            fortprintf(fd, "      end do\n");
            fortprintf(fd, "\n");
         }
         else if (dim_ptr->namelist_defined) {
            fortprintf(fd, "      read%sStart = 1\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "      read%sCount = block %% mesh %% %s\n", dim_ptr->name_in_file+1, dim_ptr->name_in_file);
            fortprintf(fd, "      allocate(send%sList)\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "      allocate(recv%sList)\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "      nullify(send%sList %% next)\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "      nullify(recv%sList %% next)\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "      recv%sList %% procID = dminfo %% my_proc_id\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "      send%sList %% procID = dminfo %% my_proc_id\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "      recv%sList %% nlist = read%sCount\n", dim_ptr->name_in_file+1, dim_ptr->name_in_file+1);
            fortprintf(fd, "      send%sList %% nlist = read%sCount\n", dim_ptr->name_in_file+1, dim_ptr->name_in_file+1);
            fortprintf(fd, "      allocate(recv%sList %% list(read%sCount))\n", dim_ptr->name_in_file+1, dim_ptr->name_in_file+1);
            fortprintf(fd, "      allocate(send%sList %% list(read%sCount))\n", dim_ptr->name_in_file+1, dim_ptr->name_in_file+1);
            fortprintf(fd, "      do i=1,read%sCount\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "         recv%sList %% list(i) = i\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "         send%sList %% list(i) = i\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "      end do\n");
            fortprintf(fd, "\n");
         }
         else {
            fortprintf(fd, "      read%sStart = 1\n", dim_ptr->name_in_code+1);
            fortprintf(fd, "      read%sCount = block %% mesh %% %s\n", dim_ptr->name_in_code+1, dim_ptr->name_in_code);
            fortprintf(fd, "      allocate(send%sList)\n", dim_ptr->name_in_code+1);
            fortprintf(fd, "      allocate(recv%sList)\n", dim_ptr->name_in_code+1);
            fortprintf(fd, "      nullify(send%sList %% next)\n", dim_ptr->name_in_code+1);
            fortprintf(fd, "      nullify(recv%sList %% next)\n", dim_ptr->name_in_code+1);
            fortprintf(fd, "      recv%sList %% procID = dminfo %% my_proc_id\n", dim_ptr->name_in_code+1);
            fortprintf(fd, "      send%sList %% procID = dminfo %% my_proc_id\n", dim_ptr->name_in_code+1);
            fortprintf(fd, "      recv%sList %% nlist = read%sCount\n", dim_ptr->name_in_code+1, dim_ptr->name_in_code+1);
            fortprintf(fd, "      send%sList %% nlist = read%sCount\n", dim_ptr->name_in_code+1, dim_ptr->name_in_code+1);
            fortprintf(fd, "      allocate(recv%sList %% list(read%sCount))\n", dim_ptr->name_in_code+1, dim_ptr->name_in_code+1);
            fortprintf(fd, "      allocate(send%sList %% list(read%sCount))\n", dim_ptr->name_in_code+1, dim_ptr->name_in_code+1);
            fortprintf(fd, "      do i=1,read%sCount\n", dim_ptr->name_in_code+1);
            fortprintf(fd, "         recv%sList %% list(i) = i\n", dim_ptr->name_in_code+1);
            fortprintf(fd, "         send%sList %% list(i) = i\n", dim_ptr->name_in_code+1);
            fortprintf(fd, "      end do\n");
            fortprintf(fd, "\n");
         }

      }

      dim_ptr = dim_ptr->next;
   }

   fclose(fd);


   /*
    *  Deallocation of exchange lists for non-decomposed fields
    */
   fd = fopen("nondecomp_dims_dealloc.inc", "w");

   dim_ptr = dims;
   while (dim_ptr) {

      if (strncmp(dim_ptr->name_in_file,"nCells",1024) != 0 &&
          strncmp(dim_ptr->name_in_file,"nEdges",1024) != 0 &&
          strncmp(dim_ptr->name_in_file,"nVertices",1024) != 0 &&
          strncmp(dim_ptr->name_in_file,"nVertLevels",1024) != 0
         ) {

         if (is_derived_dim(dim_ptr->name_in_code)) {
            fortprintf(fd, "      deallocate(recv%sList %% list)\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "      deallocate(send%sList %% list)\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "      deallocate(send%sList)\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "      deallocate(recv%sList)\n", dim_ptr->name_in_file+1);
         }
         else if (dim_ptr->constant_value > 0) {
            fortprintf(fd, "      deallocate(recv%sList %% list)\n", dim_ptr->name_in_file);
            fortprintf(fd, "      deallocate(send%sList %% list)\n", dim_ptr->name_in_file);
            fortprintf(fd, "      deallocate(send%sList)\n", dim_ptr->name_in_file);
            fortprintf(fd, "      deallocate(recv%sList)\n", dim_ptr->name_in_file);
         }
         else if (dim_ptr->namelist_defined) {
            fortprintf(fd, "      deallocate(recv%sList %% list)\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "      deallocate(send%sList %% list)\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "      deallocate(send%sList)\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "      deallocate(recv%sList)\n", dim_ptr->name_in_file+1);
         }
         else {
            fortprintf(fd, "      deallocate(recv%sList %% list)\n", dim_ptr->name_in_code+1);
            fortprintf(fd, "      deallocate(send%sList %% list)\n", dim_ptr->name_in_code+1);
            fortprintf(fd, "      deallocate(send%sList)\n", dim_ptr->name_in_code+1);
            fortprintf(fd, "      deallocate(recv%sList)\n", dim_ptr->name_in_code+1);
         }

      }

      dim_ptr = dim_ptr->next;
   }

   fclose(fd);


   /*
    *  Definitions of read bounds and exchange lists for non-decomposed fields
    */
   fd = fopen("nondecomp_outputs.inc", "w");

   dim_ptr = dims;
   while (dim_ptr) {

      if (strncmp(dim_ptr->name_in_file,"nCells",1024) != 0 &&
          strncmp(dim_ptr->name_in_file,"nEdges",1024) != 0 &&
          strncmp(dim_ptr->name_in_file,"nVertices",1024) != 0 &&
          strncmp(dim_ptr->name_in_file,"nVertLevels",1024) != 0
         ) {

         if (is_derived_dim(dim_ptr->name_in_code)) {
            fortprintf(fd, "      integer :: %sGlobal\n", dim_ptr->name_in_file);
            fortprintf(fd, "      type (exchange_list), pointer :: send%sList\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "      type (exchange_list), pointer :: recv%sList\n", dim_ptr->name_in_file+1);
         }
         else if (dim_ptr->constant_value > 0) {
            fortprintf(fd, "      type (exchange_list), pointer :: send%sList\n", dim_ptr->name_in_file);
            fortprintf(fd, "      type (exchange_list), pointer :: recv%sList\n", dim_ptr->name_in_file);
         }
         else if (dim_ptr->namelist_defined) {
            fortprintf(fd, "      integer :: %sGlobal\n", dim_ptr->name_in_file);
            fortprintf(fd, "      type (exchange_list), pointer :: send%sList\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "      type (exchange_list), pointer :: recv%sList\n", dim_ptr->name_in_file+1);
         }
         else {
            fortprintf(fd, "      integer :: %sGlobal\n", dim_ptr->name_in_code);
            fortprintf(fd, "      type (exchange_list), pointer :: send%sList\n", dim_ptr->name_in_code+1);
            fortprintf(fd, "      type (exchange_list), pointer :: recv%sList\n", dim_ptr->name_in_code+1);
         }

      }

      dim_ptr = dim_ptr->next;
   }

   fortprintf(fd, "\n");

   dim_ptr = dims;
   while (dim_ptr) {

      if (strncmp(dim_ptr->name_in_file,"nCells",1024) != 0 &&
          strncmp(dim_ptr->name_in_file,"nEdges",1024) != 0 &&
          strncmp(dim_ptr->name_in_file,"nVertices",1024) != 0 &&
          strncmp(dim_ptr->name_in_file,"nVertLevels",1024) != 0
         ) {

         if (is_derived_dim(dim_ptr->name_in_code)) {
            fortprintf(fd, "      %sGlobal = domain %% blocklist %% mesh %% %s\n", dim_ptr->name_in_file, dim_ptr->name_in_code);
            fortprintf(fd, "      if (domain %% dminfo %% my_proc_id == 0) then\n");
            fortprintf(fd, "         allocate(send%sList)\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "         allocate(recv%sList)\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "         nullify(send%sList %% next)\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "         nullify(recv%sList %% next)\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "         recv%sList %% procID = 0\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "         send%sList %% procID = 0\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "         recv%sList %% nlist = %sGlobal\n", dim_ptr->name_in_file+1, dim_ptr->name_in_file);
            fortprintf(fd, "         send%sList %% nlist = %sGlobal\n", dim_ptr->name_in_file+1, dim_ptr->name_in_file);
            fortprintf(fd, "         allocate(recv%sList %% list(%sGlobal))\n", dim_ptr->name_in_file+1, dim_ptr->name_in_file);
            fortprintf(fd, "         allocate(send%sList %% list(%sGlobal))\n", dim_ptr->name_in_file+1, dim_ptr->name_in_file);
            fortprintf(fd, "         do i=1,%sGlobal\n", dim_ptr->name_in_file);
            fortprintf(fd, "            recv%sList %% list(i) = i\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "            send%sList %% list(i) = i\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "         end do\n");
            fortprintf(fd, "      else\n");
            fortprintf(fd, "         nullify(send%sList)\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "         nullify(recv%sList)\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "      end if\n");
            fortprintf(fd, "\n");
         }
         else if (dim_ptr->constant_value > 0) {
            fortprintf(fd, "      if (domain %% dminfo %% my_proc_id == 0) then\n");
            fortprintf(fd, "         allocate(send%sList)\n", dim_ptr->name_in_file);
            fortprintf(fd, "         allocate(recv%sList)\n", dim_ptr->name_in_file);
            fortprintf(fd, "         nullify(send%sList %% next)\n", dim_ptr->name_in_file);
            fortprintf(fd, "         nullify(recv%sList %% next)\n", dim_ptr->name_in_file);
            fortprintf(fd, "         recv%sList %% procID = 0\n", dim_ptr->name_in_file);
            fortprintf(fd, "         send%sList %% procID = 0\n", dim_ptr->name_in_file);
            fortprintf(fd, "         recv%sList %% nlist = %s\n", dim_ptr->name_in_file, dim_ptr->name_in_code);
            fortprintf(fd, "         send%sList %% nlist = %s\n", dim_ptr->name_in_file, dim_ptr->name_in_code);
            fortprintf(fd, "         allocate(recv%sList %% list(%s))\n", dim_ptr->name_in_file, dim_ptr->name_in_code);
            fortprintf(fd, "         allocate(send%sList %% list(%s))\n", dim_ptr->name_in_file, dim_ptr->name_in_code);
            fortprintf(fd, "         do i=1,%s\n", dim_ptr->name_in_code);
            fortprintf(fd, "            recv%sList %% list(i) = i\n", dim_ptr->name_in_file);
            fortprintf(fd, "            send%sList %% list(i) = i\n", dim_ptr->name_in_file);
            fortprintf(fd, "         end do\n");
            fortprintf(fd, "      else\n");
            fortprintf(fd, "         nullify(send%sList)\n", dim_ptr->name_in_file);
            fortprintf(fd, "         nullify(recv%sList)\n", dim_ptr->name_in_file);
            fortprintf(fd, "      end if\n");
            fortprintf(fd, "\n");
         }
         else if (dim_ptr->namelist_defined) {
            fortprintf(fd, "      %sGlobal = domain %% blocklist %% mesh %% %s\n", dim_ptr->name_in_file, dim_ptr->name_in_file);
            fortprintf(fd, "      if (domain %% dminfo %% my_proc_id == 0) then\n");
            fortprintf(fd, "         allocate(send%sList)\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "         allocate(recv%sList)\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "         nullify(send%sList %% next)\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "         nullify(recv%sList %% next)\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "         recv%sList %% procID = 0\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "         send%sList %% procID = 0\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "         recv%sList %% nlist = %sGlobal\n", dim_ptr->name_in_file+1, dim_ptr->name_in_file);
            fortprintf(fd, "         send%sList %% nlist = %sGlobal\n", dim_ptr->name_in_file+1, dim_ptr->name_in_file);
            fortprintf(fd, "         allocate(recv%sList %% list(%sGlobal))\n", dim_ptr->name_in_file+1, dim_ptr->name_in_file);
            fortprintf(fd, "         allocate(send%sList %% list(%sGlobal))\n", dim_ptr->name_in_file+1, dim_ptr->name_in_file);
            fortprintf(fd, "         do i=1,%sGlobal\n", dim_ptr->name_in_file);
            fortprintf(fd, "            recv%sList %% list(i) = i\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "            send%sList %% list(i) = i\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "         end do\n");
            fortprintf(fd, "      else\n");
            fortprintf(fd, "         nullify(send%sList)\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "         nullify(recv%sList)\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "      end if\n");
            fortprintf(fd, "\n");
         }
         else {
            fortprintf(fd, "      %sGlobal = domain %% blocklist %% mesh %% %s\n", dim_ptr->name_in_code, dim_ptr->name_in_code);
            fortprintf(fd, "      if (domain %% dminfo %% my_proc_id == 0) then\n");
            fortprintf(fd, "         allocate(send%sList)\n", dim_ptr->name_in_code+1);
            fortprintf(fd, "         allocate(recv%sList)\n", dim_ptr->name_in_code+1);
            fortprintf(fd, "         nullify(send%sList %% next)\n", dim_ptr->name_in_code+1);
            fortprintf(fd, "         nullify(recv%sList %% next)\n", dim_ptr->name_in_code+1);
            fortprintf(fd, "         recv%sList %% procID = 0\n", dim_ptr->name_in_code+1);
            fortprintf(fd, "         send%sList %% procID = 0\n", dim_ptr->name_in_code+1);
            fortprintf(fd, "         recv%sList %% nlist = %sGlobal\n", dim_ptr->name_in_code+1, dim_ptr->name_in_code);
            fortprintf(fd, "         send%sList %% nlist = %sGlobal\n", dim_ptr->name_in_code+1, dim_ptr->name_in_code);
            fortprintf(fd, "         allocate(recv%sList %% list(%sGlobal))\n", dim_ptr->name_in_code+1, dim_ptr->name_in_code);
            fortprintf(fd, "         allocate(send%sList %% list(%sGlobal))\n", dim_ptr->name_in_code+1, dim_ptr->name_in_code);
            fortprintf(fd, "         do i=1,%sGlobal\n", dim_ptr->name_in_code);
            fortprintf(fd, "            recv%sList %% list(i) = i\n", dim_ptr->name_in_code+1);
            fortprintf(fd, "            send%sList %% list(i) = i\n", dim_ptr->name_in_code+1);
            fortprintf(fd, "         end do\n");
            fortprintf(fd, "      else\n");
            fortprintf(fd, "         nullify(send%sList)\n", dim_ptr->name_in_code+1);
            fortprintf(fd, "         nullify(recv%sList)\n", dim_ptr->name_in_code+1);
            fortprintf(fd, "      end if\n");
            fortprintf(fd, "\n");
         }

      }

      dim_ptr = dim_ptr->next;
   }

   fclose(fd);


   /*
    *  Deallocation of exchange lists for non-decomposed fields
    */
   fd = fopen("nondecomp_outputs_dealloc.inc", "w");

   dim_ptr = dims;
   while (dim_ptr) {

      if (strncmp(dim_ptr->name_in_file,"nCells",1024) != 0 &&
          strncmp(dim_ptr->name_in_file,"nEdges",1024) != 0 &&
          strncmp(dim_ptr->name_in_file,"nVertices",1024) != 0 &&
          strncmp(dim_ptr->name_in_file,"nVertLevels",1024) != 0
         ) {

         if (is_derived_dim(dim_ptr->name_in_code)) {
            fortprintf(fd, "      if (domain %% dminfo %% my_proc_id == 0) then\n");
            fortprintf(fd, "         deallocate(recv%sList %% list)\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "         deallocate(send%sList %% list)\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "         deallocate(send%sList)\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "         deallocate(recv%sList)\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "      end if\n");
         }
         else if (dim_ptr->constant_value > 0) {
            fortprintf(fd, "      if (domain %% dminfo %% my_proc_id == 0) then\n");
            fortprintf(fd, "         deallocate(recv%sList %% list)\n", dim_ptr->name_in_file);
            fortprintf(fd, "         deallocate(send%sList %% list)\n", dim_ptr->name_in_file);
            fortprintf(fd, "         deallocate(send%sList)\n", dim_ptr->name_in_file);
            fortprintf(fd, "         deallocate(recv%sList)\n", dim_ptr->name_in_file);
            fortprintf(fd, "      end if\n");
         }
         else if (dim_ptr->namelist_defined) {
            fortprintf(fd, "      if (domain %% dminfo %% my_proc_id == 0) then\n");
            fortprintf(fd, "         deallocate(recv%sList %% list)\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "         deallocate(send%sList %% list)\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "         deallocate(send%sList)\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "         deallocate(recv%sList)\n", dim_ptr->name_in_file+1);
            fortprintf(fd, "      end if\n");
         }
         else {
            fortprintf(fd, "      if (domain %% dminfo %% my_proc_id == 0) then\n");
            fortprintf(fd, "         deallocate(recv%sList %% list)\n", dim_ptr->name_in_code+1);
            fortprintf(fd, "         deallocate(send%sList %% list)\n", dim_ptr->name_in_code+1);
            fortprintf(fd, "         deallocate(send%sList)\n", dim_ptr->name_in_code+1);
            fortprintf(fd, "         deallocate(recv%sList)\n", dim_ptr->name_in_code+1);
            fortprintf(fd, "      end if\n");
         }

      }

      dim_ptr = dim_ptr->next;
   }

   fclose(fd);
   

   /*
    *  Generate read and distribute code
    */
   fd = fopen("io_input_fields.inc", "w");

   group_ptr = groups;
   while (group_ptr) {
      var_list_ptr = group_ptr->vlist;
      while (var_list_ptr) {
         var_ptr = var_list_ptr->var;

         if (group_ptr->vlist->var->ntime_levs > 1)
            snprintf(struct_deref, 1024, "block %% %s %% time_levs(1) %% %s", group_ptr->name, group_ptr->name);
         else
            snprintf(struct_deref, 1024, "block %% %s", group_ptr->name);

         i = 1;
         dimlist_ptr = var_ptr->dimlist;
         if (var_ptr->vtype == INTEGER) sprintf(vtype, "int"); 
         else if (var_ptr->vtype == REAL) sprintf(vtype, "real"); 
         else if (var_ptr->vtype == CHARACTER) sprintf(vtype, "char"); 
   
         if (strncmp(var_ptr->super_array, "-", 1024) != 0) {
            fortprintf(fd, "      if ((%s %% %s %% ioinfo %% input .and. input_obj %% stream == STREAM_INPUT) .or. &\n", struct_deref, var_ptr->super_array);
            fortprintf(fd, "          (%s %% %s %% ioinfo %% sfc .and. input_obj %% stream == STREAM_SFC) .or. &\n", struct_deref, var_ptr->super_array);
            fortprintf(fd, "          (%s %% %s %% ioinfo %% restart .and. input_obj %% stream == STREAM_RESTART)) then\n", struct_deref, var_ptr->super_array);
         }
         else {
            fortprintf(fd, "      if ((%s %% %s %% ioinfo %% input .and. input_obj %% stream == STREAM_INPUT) .or. &\n", struct_deref, var_ptr->name_in_code);
            fortprintf(fd, "          (%s %% %s %% ioinfo %% sfc .and. input_obj %% stream == STREAM_SFC) .or. &\n", struct_deref, var_ptr->name_in_code);
            fortprintf(fd, "          (%s %% %s %% ioinfo %% restart .and. input_obj %% stream == STREAM_RESTART)) then\n", struct_deref, var_ptr->name_in_code);
         }
         while (dimlist_ptr) {
               if (i < var_ptr->ndims) {
                  fortprintf(fd, "      %s%id %% ioinfo %% start(%i) = 1\n", vtype, var_ptr->ndims, i);
                  if (dimlist_ptr->dim->constant_value < 0)
                     if (!dimlist_ptr->dim->namelist_defined) fortprintf(fd, "      %s%id %% ioinfo %% count(%i) = block %% mesh %% %s\n", vtype, var_ptr->ndims, i, dimlist_ptr->dim->name_in_code);
                     else fortprintf(fd, "      %s%id %% ioinfo %% count(%i) = block %% mesh %% %s\n", vtype, var_ptr->ndims, i, dimlist_ptr->dim->name_in_file);
                  else
                     fortprintf(fd, "      %s%id %% ioinfo %% count(%i) = %s\n", vtype, var_ptr->ndims, i, dimlist_ptr->dim->name_in_code);
               }
               else {
                  if (is_derived_dim(dimlist_ptr->dim->name_in_code)) {
                     fortprintf(fd, "      %s%id %% ioinfo %% start(%i) = read%sStart\n", vtype, var_ptr->ndims, i, dimlist_ptr->dim->name_in_file+1);
                     fortprintf(fd, "      %s%id %% ioinfo %% count(%i) = read%sCount\n", vtype, var_ptr->ndims, i, dimlist_ptr->dim->name_in_file+1);
                  }
                  else if (dimlist_ptr->dim->constant_value > 0) {
                     fortprintf(fd, "      %s%id %% ioinfo %% start(%i) = read%sStart\n", vtype, var_ptr->ndims, i, dimlist_ptr->dim->name_in_file);
                     fortprintf(fd, "      %s%id %% ioinfo %% count(%i) = read%sCount\n", vtype, var_ptr->ndims, i, dimlist_ptr->dim->name_in_file);
                  }
                  else {
                     if (dimlist_ptr->dim->namelist_defined) {
                        fortprintf(fd, "      %s%id %% ioinfo %% start(%i) = read%sStart\n", vtype, var_ptr->ndims, i, dimlist_ptr->dim->name_in_file+1);
                        fortprintf(fd, "      %s%id %% ioinfo %% count(%i) = read%sCount\n", vtype, var_ptr->ndims, i, dimlist_ptr->dim->name_in_file+1);
                     }
                     else {
                        fortprintf(fd, "      %s%id %% ioinfo %% start(%i) = read%sStart\n", vtype, var_ptr->ndims, i, dimlist_ptr->dim->name_in_code+1);
                        fortprintf(fd, "      %s%id %% ioinfo %% count(%i) = read%sCount\n", vtype, var_ptr->ndims, i, dimlist_ptr->dim->name_in_code+1);
                     }
                  }
               }
            dimlist_ptr = dimlist_ptr->next;
            i++;
         }
   
         if (var_ptr->ndims > 0) {
            fortprintf(fd, "      allocate(%s%id %% array(", vtype, var_ptr->ndims);
            i = 1;
            dimlist_ptr = var_ptr->dimlist;
      
            if (i < var_ptr->ndims) {
               if (dimlist_ptr->dim->constant_value < 0)
                  if (!dimlist_ptr->dim->namelist_defined) fortprintf(fd, "block %% mesh %% %s", dimlist_ptr->dim->name_in_code);
                  else fortprintf(fd, "block %% mesh %% %s", dimlist_ptr->dim->name_in_file);
               else
                  fortprintf(fd, "%s", dimlist_ptr->dim->name_in_code);
            }
            else {
               if (is_derived_dim(dimlist_ptr->dim->name_in_code)) 
                  fortprintf(fd, "read%sCount", dimlist_ptr->dim->name_in_file+1);
               else if (dimlist_ptr->dim->constant_value > 0) 
                  fortprintf(fd, "read%sCount", dimlist_ptr->dim->name_in_file);
               else
                  if (dimlist_ptr->dim->namelist_defined) fortprintf(fd, "read%sCount", dimlist_ptr->dim->name_in_file+1);
                  else fortprintf(fd, "read%sCount", dimlist_ptr->dim->name_in_code+1);
            }
       
            dimlist_ptr = dimlist_ptr->next;
            i++;
            while (dimlist_ptr) {
               if (i < var_ptr->ndims) {
                  if (dimlist_ptr->dim->constant_value < 0)
                     if (!dimlist_ptr->dim->namelist_defined) fortprintf(fd, ", block %% mesh %% %s", dimlist_ptr->dim->name_in_code);
                     else fortprintf(fd, ", block %% mesh %% %s", dimlist_ptr->dim->name_in_file);
                  else
                     fortprintf(fd, ", %s", dimlist_ptr->dim->name_in_code);
               }
               else {
                  if (is_derived_dim(dimlist_ptr->dim->name_in_code))
                     fortprintf(fd, ", read%sCount", dimlist_ptr->dim->name_in_file+1);
                  else if (dimlist_ptr->dim->constant_value > 0)
                     fortprintf(fd, ", read%sCount", dimlist_ptr->dim->name_in_file);
                  else
                     if (!dimlist_ptr->dim->namelist_defined) fortprintf(fd, ", read%sCount", dimlist_ptr->dim->name_in_code+1);
                     else fortprintf(fd, ", read%sCount", dimlist_ptr->dim->name_in_file+1);
               }
               dimlist_ptr = dimlist_ptr->next;
               i++;
            }
            fortprintf(fd, "))\n\n");
   
            if (strncmp(var_ptr->super_array, "-", 1024) != 0) {
               fortprintf(fd, "      allocate(super_%s%id(", vtype, var_ptr->ndims);
               i = 1;
               dimlist_ptr = var_ptr->dimlist;
         
               if (i < var_ptr->ndims) {
                  if (dimlist_ptr->dim->constant_value < 0)
                     if (!dimlist_ptr->dim->namelist_defined) fortprintf(fd, "block %% mesh %% %s", dimlist_ptr->dim->name_in_code);
                     else fortprintf(fd, "block %% mesh %% %s", dimlist_ptr->dim->name_in_file);
                  else
                     if (!dimlist_ptr->dim->namelist_defined) fortprintf(fd, "%s", dimlist_ptr->dim->name_in_code);
                     else fortprintf(fd, "%s", dimlist_ptr->dim->name_in_file);
               }
               dimlist_ptr = dimlist_ptr->next;
               i++;
               while (dimlist_ptr) {
                  if (dimlist_ptr->dim->constant_value < 0)
                     if (!dimlist_ptr->dim->namelist_defined) fortprintf(fd, ", block %% mesh %% %s", dimlist_ptr->dim->name_in_code);
                     else fortprintf(fd, ", block %% mesh %% %s", dimlist_ptr->dim->name_in_file);
                  else
                     if (!dimlist_ptr->dim->namelist_defined) fortprintf(fd, ", %s", dimlist_ptr->dim->name_in_code);
                     else fortprintf(fd, ", %s", dimlist_ptr->dim->name_in_file);

                  dimlist_ptr = dimlist_ptr->next;
                  i++;
               }
               fortprintf(fd, "))\n\n");
            }
         }
   
         fortprintf(fd, "      %s%id %% ioinfo %% fieldName = \'%s\'\n", vtype, var_ptr->ndims, var_ptr->name_in_file);
         if (var_ptr->timedim)
            fortprintf(fd, "      call mpas_io_input_field_time(input_obj, %s%id)\n", vtype, var_ptr->ndims);
         else
            fortprintf(fd, "      call mpas_io_input_field(input_obj, %s%id)\n", vtype, var_ptr->ndims);
   
         if (var_ptr->ndims > 0) {
            fortprintf(fd, "      call mpas_dmpar_alltoall_field(dminfo, &\n");
            if (strncmp(var_ptr->super_array, "-", 1024) != 0)
               fortprintf(fd, "                                %s%id %% array, super_%s%id, &\n", vtype, var_ptr->ndims, vtype, var_ptr->ndims);
            else
               fortprintf(fd, "                                %s%id %% array, %s %% %s %% array, &\n", vtype, var_ptr->ndims, struct_deref, var_ptr->name_in_code);
      
            i = 1;
            dimlist_ptr = var_ptr->dimlist;
            
            if (i < var_ptr->ndims)
               if (dimlist_ptr->dim->constant_value < 0)
                  if (!dimlist_ptr->dim->namelist_defined) fortprintf(fd, "                                block %% mesh %% %s", dimlist_ptr->dim->name_in_code);
                  else fortprintf(fd, "                                block %% mesh %% %s", dimlist_ptr->dim->name_in_file);
               else
                  fortprintf(fd, "                                %s", dimlist_ptr->dim->name_in_code);
            else {
               lastdim = dimlist_ptr;
               if (is_derived_dim(dimlist_ptr->dim->name_in_code)) {
                  fortprintf(fd, "                                read%sCount", dimlist_ptr->dim->name_in_file+1);
               }
               else if (dimlist_ptr->dim->constant_value > 0) {
                  fortprintf(fd, "                                read%sCount", dimlist_ptr->dim->name_in_file);
               }
               else
                  if (!dimlist_ptr->dim->namelist_defined) fortprintf(fd, "                                read%sCount", dimlist_ptr->dim->name_in_code+1);
                  else fortprintf(fd, "                                read%sCount", dimlist_ptr->dim->name_in_file+1);
            }
       
            dimlist_ptr = dimlist_ptr->next;
            i++;
            while (dimlist_ptr) {
               if (i < var_ptr->ndims)
                  if (dimlist_ptr->dim->constant_value < 0)
                     if (!dimlist_ptr->dim->namelist_defined) fortprintf(fd, ", block %% mesh %% %s", dimlist_ptr->dim->name_in_code);
                     else fortprintf(fd, ", block %% mesh %% %s", dimlist_ptr->dim->name_in_file);
                  else
                     fortprintf(fd, ", %s", dimlist_ptr->dim->name_in_code);
               else {
                  lastdim = dimlist_ptr;
                  if (is_derived_dim(dimlist_ptr->dim->name_in_code)) {
                     fortprintf(fd, ", read%sCount", dimlist_ptr->dim->name_in_file+1);
                  }
                  else if (dimlist_ptr->dim->constant_value > 0) {
                     fortprintf(fd, ", read%sCount", dimlist_ptr->dim->name_in_file);
                  }
                  else
                     if (!dimlist_ptr->dim->namelist_defined) fortprintf(fd, ", read%sCount", dimlist_ptr->dim->name_in_code+1);
                     else fortprintf(fd, ", read%sCount", dimlist_ptr->dim->name_in_file+1);
               }
               dimlist_ptr = dimlist_ptr->next;
               i++;
            }
            if (lastdim->dim->namelist_defined) fortprintf(fd, ", block %% mesh %% %s, &\n", lastdim->dim->name_in_file);
            else if (lastdim->dim->constant_value > 0) fortprintf(fd, ", %s, &\n", lastdim->dim->name_in_code);
            else fortprintf(fd, ", block %% mesh %% %s, &\n", lastdim->dim->name_in_code);
      
            if (is_derived_dim(lastdim->dim->name_in_code)) 
               fortprintf(fd, "                                send%sList, recv%sList)\n", lastdim->dim->name_in_file+1, lastdim->dim->name_in_file+1);
            else if (lastdim->dim->constant_value > 0) 
               fortprintf(fd, "                                send%sList, recv%sList)\n", lastdim->dim->name_in_file, lastdim->dim->name_in_file);
            else
               if (lastdim->dim->namelist_defined) 
                  fortprintf(fd, "                                send%sList, recv%sList)\n", lastdim->dim->name_in_file+1, lastdim->dim->name_in_file+1);
               else
                  fortprintf(fd, "                                send%sList, recv%sList)\n", lastdim->dim->name_in_code+1, lastdim->dim->name_in_code+1);
   
   
            /* Copy from super_ array to field */
            if (strncmp(var_ptr->super_array, "-", 1024) != 0) {
               i = 1;
               dimlist_ptr = var_ptr->dimlist;
               while (i <= var_ptr->ndims) {
                  if (dimlist_ptr->dim->constant_value < 0)
                     if (!dimlist_ptr->dim->namelist_defined) fortprintf(fd, "      do i%i=1,block %% mesh %% %s\n", i, dimlist_ptr->dim->name_in_code);
                     else fortprintf(fd, "      do i%i=1,block %% mesh %% %s\n", i, dimlist_ptr->dim->name_in_file);
                  else
                     fortprintf(fd, "      do i%i=1,%s\n", i, dimlist_ptr->dim->name_in_code);
      
                  i++;
                  dimlist_ptr = dimlist_ptr->next;
               }
      
               fortprintf(fd, "         %s %% %s %% array(%s %% index_%s,", struct_deref, var_ptr->super_array, struct_deref, var_ptr->name_in_code);

               for(i=1; i<=var_ptr->ndims; i++) {
                  fortprintf(fd, "i%i",i);
                  if (i < var_ptr->ndims) fortprintf(fd, ",");
               }
               fortprintf(fd, ") = super_%s%id(", vtype, var_ptr->ndims);
               for(i=1; i<=var_ptr->ndims; i++) {
                  fortprintf(fd, "i%i",i);
                  if (i < var_ptr->ndims) fortprintf(fd, ",");
               }
               fortprintf(fd, ")\n");
      
               i = 1;
               while (i <= var_ptr->ndims) {
                  fortprintf(fd, "      end do\n");
                  i++;
               }
            }
   
            fortprintf(fd, "      deallocate(%s%id %% array)\n", vtype, var_ptr->ndims);
            if (strncmp(var_ptr->super_array, "-", 1024) != 0)
               fortprintf(fd, "      deallocate(super_%s%id)\n", vtype, var_ptr->ndims);
         }
         else {
            fortprintf(fd, "      %s %% %s %% scalar = %s%id %% scalar\n", struct_deref, var_ptr->name_in_code, vtype, var_ptr->ndims);
         }
        
         fortprintf(fd, "      end if\n\n");
   
         var_list_ptr = var_list_ptr->next;
      }
      group_ptr = group_ptr->next;
   }

   fclose(fd);
#endif


   /*
    * MGD NEW CODE
    */
   fd = fopen("add_input_fields.inc", "w");

   group_ptr = groups;
   while (group_ptr) {
      var_list_ptr = group_ptr->vlist;
      while (var_list_ptr) {
         var_ptr = var_list_ptr->var;

         if (var_ptr->ntime_levs > 1)
            snprintf(struct_deref, 1024, "blocklist %% %s %% time_levs(1) %% %s", group_ptr->name, group_ptr->name);
         else
            snprintf(struct_deref, 1024, "blocklist %% %s", group_ptr->name);
         
         if (strncmp(var_ptr->super_array, "-", 1024) != 0) {
            fortprintf(fd, "      if ((%s %% %s %% ioinfo %% input .and. input_obj %% stream == STREAM_INPUT) .or. &\n", struct_deref, var_ptr->super_array);
            fortprintf(fd, "          (%s %% %s %% ioinfo %% restart .and. input_obj %% stream == STREAM_RESTART) .or. &\n", struct_deref, var_ptr->super_array);
            fortprintf(fd, "          (%s %% %s %% ioinfo %% sfc .and. input_obj %% stream == STREAM_SFC)) then\n", struct_deref, var_ptr->super_array);
            memcpy(super_array, var_ptr->super_array, 1024);
/*            fortprintf(fd, "         write(0,*) \'adding input field %s\'\n", var_ptr->super_array); */
            fortprintf(fd, "         call MPAS_streamAddField(input_obj %% io_stream, %s %% %s, nferr)\n", struct_deref, var_ptr->super_array);
            while (var_list_ptr && strncmp(super_array, var_list_ptr->var->super_array, 1024) == 0) {
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


   /*
    * MGD NEW CODE
    */
   fd = fopen("exchange_input_field_halos.inc", "w");
   fd2 = fopen("non_decomp_copy_input_fields.inc", "w");

   group_ptr = groups;
   while (group_ptr) {
      var_list_ptr = group_ptr->vlist;
      while (var_list_ptr) {
         var_ptr = var_list_ptr->var;

         dimlist_ptr = var_ptr->dimlist;
         i = 1;
		 if(var_ptr->persistence == PERSISTENT){
         while (dimlist_ptr) {
            if (i == var_ptr->ndims) { 

                  if (var_ptr->ntime_levs > 1) {
                     snprintf(struct_deref, 1024, "domain %% blocklist %% %s %% time_levs(1) %% %s", group_ptr->name, group_ptr->name);
				  } else {
                     snprintf(struct_deref, 1024, "domain %% blocklist %% %s", group_ptr->name);
				  }

               if (!strncmp(dimlist_ptr->dim->name_in_file, "nCells", 1024) ||
                   !strncmp(dimlist_ptr->dim->name_in_file, "nEdges", 1024) ||
                   !strncmp(dimlist_ptr->dim->name_in_file, "nVertices", 1024)) {
                  
                  if (strncmp(var_ptr->super_array, "-", 1024) != 0) {
                     fortprintf(fd, "      if ((%s %% %s %% ioinfo %% input .and. input_obj %% stream == STREAM_INPUT) .or. &\n", struct_deref, var_ptr->super_array);
                     fortprintf(fd, "          (%s %% %s %% ioinfo %% restart .and. input_obj %% stream == STREAM_RESTART) .or. &\n", struct_deref, var_ptr->super_array);
                     fortprintf(fd, "          (%s %% %s %% ioinfo %% sfc .and. input_obj %% stream == STREAM_SFC)) then\n", struct_deref, var_ptr->super_array);
                     memcpy(super_array, var_ptr->super_array, 1024);
/*                     fortprintf(fd, "         write(0,*) \'exchange halo for %s\'\n", var_ptr->super_array); */
                     fortprintf(fd, "         call mpas_dmpar_exch_halo_field(%s %% %s)\n", struct_deref, var_ptr->super_array);
                     while (var_list_ptr && strncmp(super_array, var_list_ptr->var->super_array, 1024) == 0) {
						var_list_ptr2 = var_list_ptr;
                        var_list_ptr = var_list_ptr->next;
                     }
					 var_list_ptr = var_list_ptr2;
                  }
                  else {
                     fortprintf(fd, "      if ((%s %% %s %% ioinfo %% input .and. input_obj %% stream == STREAM_INPUT) .or. &\n", struct_deref, var_ptr->name_in_code);
                     fortprintf(fd, "          (%s %% %s %% ioinfo %% restart .and. input_obj %% stream == STREAM_RESTART) .or. &\n", struct_deref, var_ptr->name_in_code);
                     fortprintf(fd, "          (%s %% %s %% ioinfo %% sfc .and. input_obj %% stream == STREAM_SFC)) then\n", struct_deref, var_ptr->name_in_code);
/*                     fortprintf(fd, "         write(0,*) \'exchange halo for %s\'\n", var_ptr->name_in_code); */
                     fortprintf(fd, "         call mpas_dmpar_exch_halo_field(%s %% %s)\n", struct_deref, var_ptr->name_in_code);
                  }
            
                  fortprintf(fd, "      end if\n\n");
   
               } else {
                  fortprintf(fd2, "      if ((%s %% %s %% ioinfo %% input .and. input_obj %% stream == STREAM_INPUT) .or. &\n", struct_deref, var_ptr->name_in_code);
                  fortprintf(fd2, "          (%s %% %s %% ioinfo %% restart .and. input_obj %% stream == STREAM_RESTART) .or. &\n", struct_deref, var_ptr->name_in_code);
                  fortprintf(fd2, "          (%s %% %s %% ioinfo %% sfc .and. input_obj %% stream == STREAM_SFC)) then\n", struct_deref, var_ptr->name_in_code);
				  fortprintf(fd2, "          call mpas_dmpar_copy_field(%s %% %s)\n", struct_deref, var_ptr->name_in_code);
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


#ifdef LEGACY_CODE
   /*
    *  Generate NetCDF reads of dimension and variable IDs
    */
   fd = fopen("netcdf_read_ids.inc", "w");

   fortprintf(fd, "      nferr = nf_inq_unlimdim(input_obj %% rd_ncid, input_obj %% rdDimIDTime)\n");
   fortprintf(fd, "      nferr = nf_inq_dimlen(input_obj %% rd_ncid, input_obj %% rdDimIDTime, input_obj %% rdLocalTime)\n");
   dim_ptr = dims;
   while (dim_ptr) {
      if (dim_ptr->constant_value < 0 && !dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) {
         fortprintf(fd, "      nferr = nf_inq_dimid(input_obj %% rd_ncid, \'%s\', input_obj %% rdDimID%s)\n", dim_ptr->name_in_file, dim_ptr->name_in_file);
         fortprintf(fd, "      nferr = nf_inq_dimlen(input_obj %% rd_ncid, input_obj %% rdDimID%s, input_obj %% rdLocal%s)\n", dim_ptr->name_in_file, dim_ptr->name_in_file);
      }
      dim_ptr = dim_ptr->next;
   }
   fortprintf(fd, "\n");

   var_ptr = vars;
   while (var_ptr) {
      fortprintf(fd, "      nferr = nf_inq_varid(input_obj %% rd_ncid, \'%s\', input_obj %% rdVarID%s)\n", var_ptr->name_in_file, var_ptr->name_in_file);
      var_ptr = var_ptr->next;
   }

   fclose(fd);


   /*
    *  Generate code to return dimension given its name
    */
   fd = fopen("get_dimension_by_name.inc", "w");

   dim_ptr = dims;
   while (dim_ptr->constant_value >= 0 || is_derived_dim(dim_ptr->name_in_code)) dim_ptr = dim_ptr->next;
   if (!dim_ptr->namelist_defined) {
      fortprintf(fd, "      if (trim(dimname) == \'%s\') then\n", dim_ptr->name_in_code);
      fortprintf(fd, "         dimsize = input_obj %% rdLocal%s\n", dim_ptr->name_in_file);
   }
   else {
      fortprintf(fd, "      if (trim(dimname) == \'%s\') then\n", dim_ptr->name_in_file);
      fortprintf(fd, "         dimsize = %s\n", dim_ptr->name_in_code);
   }
   dim_ptr = dim_ptr->next;
   while (dim_ptr) {
      if (dim_ptr->constant_value < 0 && !is_derived_dim(dim_ptr->name_in_code)) {
         if (!dim_ptr->namelist_defined) {
            fortprintf(fd, "      else if (trim(dimname) == \'%s\') then\n", dim_ptr->name_in_code);
            fortprintf(fd, "         dimsize = input_obj %% rdLocal%s\n", dim_ptr->name_in_file);
         }
         else {
            fortprintf(fd, "      else if (trim(dimname) == \'%s\') then\n", dim_ptr->name_in_file);
            fortprintf(fd, "         dimsize = %s\n", dim_ptr->name_in_code);
         }
      }
      dim_ptr = dim_ptr->next;
   }
   fortprintf(fd, "      end if\n");

   fclose(fd);
   
   
   /*
    *  Generate code to read 0d, 1d, 2d, 3d time-invariant fields
    */
   for(j=0; j<3; j++) {
      for(i=0; i<=3; i++) {
         if (j == 0) {
            sprintf(fname, "input_field%idinteger.inc", i);
            ivtype = INTEGER;
         }
         else if (j == 1) {
            sprintf(fname, "input_field%idreal.inc", i);
            ivtype = REAL;
         }
         else if (j == 2) {
            sprintf(fname, "input_field%idchar.inc", i);
            ivtype = CHARACTER;
         }
         fd = fopen(fname, "w");
   
         var_ptr = vars;
         while (var_ptr && (var_ptr->ndims != i || var_ptr->vtype != ivtype || var_ptr->timedim)) var_ptr = var_ptr->next;
         if (var_ptr) {
            fortprintf(fd, "      if (trim(field %% ioinfo %% fieldName) == \'%s\') then\n", var_ptr->name_in_file);
            fortprintf(fd, "         varID = input_obj %% rdVarID%s\n", var_ptr->name_in_file);
            var_ptr = var_ptr->next;
            while (var_ptr) {
               if (var_ptr->ndims == i && var_ptr->vtype == ivtype && !var_ptr->timedim) {
                  fortprintf(fd, "      else if (trim(field %% ioinfo %% fieldName) == \'%s\') then\n", var_ptr->name_in_file);
                  fortprintf(fd, "         varID = input_obj %% rdVarID%s\n", var_ptr->name_in_file);
               }
               var_ptr = var_ptr->next;
            }
            fortprintf(fd, "      end if\n");
         }
      
         fclose(fd);
      } 
   } 
   
   
   /*
    *  Generate code to read 0d, 1d, 2d, 3d time-varying fields
    */
   for(j=0; j<2; j++) {
      for(i=0; i<=3; i++) { 
         if (j == 0) {
            sprintf(fname, "input_field%idinteger_time.inc", i);
            ivtype = INTEGER;
         }
         else {
            sprintf(fname, "input_field%idreal_time.inc", i);
            ivtype = REAL;
         }
         fd = fopen(fname, "w");
      
         var_ptr = vars;
         while (var_ptr && (var_ptr->ndims != i || var_ptr->vtype != ivtype || !var_ptr->timedim)) var_ptr = var_ptr->next;
         if (var_ptr) {
            fortprintf(fd, "      if (trim(field %% ioinfo %% fieldName) == \'%s\') then\n", var_ptr->name_in_file);
            fortprintf(fd, "         varID = input_obj %% rdVarID%s\n", var_ptr->name_in_file);
            var_ptr = var_ptr->next;
            while (var_ptr) {
               if (var_ptr->ndims == i && var_ptr->vtype == ivtype && var_ptr->timedim) {
                  fortprintf(fd, "      else if (trim(field %% ioinfo %% fieldName) == \'%s\') then\n", var_ptr->name_in_file);
                  fortprintf(fd, "         varID = input_obj %% rdVarID%s\n", var_ptr->name_in_file);
               }
               var_ptr = var_ptr->next;
   
	    }
	    fortprintf(fd, "      end if\n");
	 }
   
	 fclose(fd);
      }
   } 
   
   
   /*
    *  Generate code to read 0d and 1d time-varying character fields
    */
   for(i=0; i<=1; i++) { 
      sprintf(fname, "input_field%idchar_time.inc", i);
      fd = fopen(fname, "w");
   
      var_ptr = vars;
      while (var_ptr && (var_ptr->ndims != i || var_ptr->vtype != CHARACTER || !var_ptr->timedim)) var_ptr = var_ptr->next;
      if (var_ptr) {
         fortprintf(fd, "      if (trim(field %% ioinfo %% fieldName) == \'%s\') then\n", var_ptr->name_in_file);
         fortprintf(fd, "         varID = input_obj %% rdVarID%s\n", var_ptr->name_in_file);
         var_ptr = var_ptr->next;
         while (var_ptr) {
            if (var_ptr->ndims == i && var_ptr->vtype == CHARACTER && var_ptr->timedim) {
               fortprintf(fd, "      else if (trim(field %% ioinfo %% fieldName) == \'%s\') then\n", var_ptr->name_in_file);
               fortprintf(fd, "         varID = input_obj %% rdVarID%s\n", var_ptr->name_in_file);
            }
            var_ptr = var_ptr->next;
         }
         fortprintf(fd, "      end if\n");
      }
   
      fclose(fd);
   } 
#endif
   
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
   char super_array[1024];
   char * cp1, * cp2;
   int i, j;
   int ivtype;
   
   
#ifdef LEGACY_CODE
   /*
    *  Generate declarations of IDs belonging in io_output_object
    */
   fd = fopen("io_output_obj_decls.inc", "w");

   fortprintf(fd, "      integer :: wrDimIDTime\n");
   dim_ptr = dims;
   while (dim_ptr) {
      fortprintf(fd, "      integer :: wrDimID%s\n", dim_ptr->name_in_file);
      dim_ptr = dim_ptr->next;
   }
   fortprintf(fd, "\n");

   var_ptr = vars;
   while (var_ptr) {
      fortprintf(fd, "      integer :: wrVarID%s\n", var_ptr->name_in_file);
      var_ptr = var_ptr->next;
   }

   fclose(fd);


   /*
    *  Generate declarations of temporary dimension variables used for arguments
    */
   fd = fopen("output_dim_actual_decls.inc", "w");

   dim_ptr = dims;
   while (dim_ptr) {
      if (dim_ptr->constant_value < 0 && !dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) fortprintf(fd, "      integer :: %sGlobal\n", dim_ptr->name_in_code);
      if (dim_ptr->constant_value < 0 && dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) fortprintf(fd, "      integer :: %sGlobal\n", dim_ptr->name_in_file);
      dim_ptr = dim_ptr->next;
   }
   fortprintf(fd, "\n");

   fclose(fd);


   /*
    *  Generate initialization of temporary dimension variables used for arguments
    */
   fd = fopen("output_dim_inits.inc", "w");

   dim_ptr = dims;
   while (dim_ptr) {
      if (dim_ptr->constant_value < 0 && !dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) fortprintf(fd, "      %sGlobal = block_ptr %% mesh %% %s\n", dim_ptr->name_in_code, dim_ptr->name_in_code);
      if (dim_ptr->constant_value < 0 && dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) fortprintf(fd, "      %sGlobal = block_ptr %% mesh %% %s\n", dim_ptr->name_in_file, dim_ptr->name_in_file);
      dim_ptr = dim_ptr->next;
   }
   fortprintf(fd, "\n");

   fclose(fd);


   /*
    *  Generate actual dimension argument list
    */
   fd = fopen("output_dim_actual_args.inc", "w");
   dim_ptr = dims;
   if (dim_ptr && dim_ptr->constant_value < 0 && !is_derived_dim(dim_ptr->name_in_code)) {
      if (!dim_ptr->namelist_defined) fortprintf(fd, "                            %sGlobal", dim_ptr->name_in_code);
      else fortprintf(fd, "                            %sGlobal", dim_ptr->name_in_file);
      dim_ptr = dim_ptr->next;
   }
   while (dim_ptr) {
      if (dim_ptr->constant_value < 0 && !dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) fortprintf(fd, ", %sGlobal", dim_ptr->name_in_code);
      if (dim_ptr->constant_value < 0 && dim_ptr->namelist_defined && !is_derived_dim(dim_ptr->name_in_code)) fortprintf(fd, ", %sGlobal", dim_ptr->name_in_file);
      dim_ptr = dim_ptr->next;
   }
   fortprintf(fd, " &\n");

   fclose(fd);


   /*
    *  Generate NetCDF calls to define dimensions, variables, and global attributes
    */
   fd = fopen("netcdf_def_dims_vars.inc", "w");

   fortprintf(fd, "      nferr = nf_def_dim(output_obj %% wr_ncid, \'Time\', NF_UNLIMITED, output_obj %% wrDimIDTime)\n");
   dim_ptr = dims;
   while (dim_ptr) {
      fortprintf(fd, "      nferr = nf_def_dim(output_obj %% wr_ncid, \'%s\', %s, output_obj %% wrDimID%s)\n", dim_ptr->name_in_file, dim_ptr->name_in_code, dim_ptr->name_in_file);
      dim_ptr = dim_ptr->next;
   }
   fortprintf(fd, "\n");

   var_ptr = vars;
   while (var_ptr) {
      fortprintf(fd, "      if (.false. &\n");
      if (var_ptr->iostreams & RESTART0) fortprintf(fd, "          .or. output_obj %% stream == RESTART &\n");
      if (var_ptr->iostreams & OUTPUT0)  fortprintf(fd, "          .or. output_obj %% stream == OUTPUT &\n");
      if (var_ptr->iostreams & SFC0)     fortprintf(fd, "          .or. output_obj %% stream == SFC &\n");
      fortprintf(fd, "      ) then\n");
      dimlist_ptr = var_ptr->dimlist;
      i = 1;
      if (var_ptr->vtype == CHARACTER)
         fortprintf(fd, "      dimlist(%i) = output_obj %% wrDimIDStrLen\n", i++);
      while(dimlist_ptr) {
         fortprintf(fd, "      dimlist(%i) = output_obj %% wrDimID%s\n", i++, dimlist_ptr->dim->name_in_file);
         dimlist_ptr = dimlist_ptr->next;
      }
      if (var_ptr->timedim) fortprintf(fd, "      dimlist(%i) = output_obj %% wrDimIDTime\n", i++);
      if (var_ptr->vtype == INTEGER)
         fortprintf(fd, "      nferr = nf_def_var(output_obj %% wr_ncid, \'%s\', NF_INT, %i, dimlist, output_obj %% wrVarID%s)\n", var_ptr->name_in_file, var_ptr->ndims + var_ptr->timedim, var_ptr->name_in_file);
      else if (var_ptr->vtype == REAL)
         fortprintf(fd, "      nferr = nf_def_var(output_obj %% wr_ncid, \'%s\', NF_DOUBLE, %i, dimlist, output_obj %% wrVarID%s)\n", var_ptr->name_in_file, var_ptr->ndims + var_ptr->timedim, var_ptr->name_in_file);
      else if (var_ptr->vtype == CHARACTER)
         fortprintf(fd, "      nferr = nf_def_var(output_obj %% wr_ncid, \'%s\', NF_CHAR, %i, dimlist, output_obj %% wrVarID%s)\n", var_ptr->name_in_file, var_ptr->ndims + var_ptr->timedim + 1, var_ptr->name_in_file);

      fortprintf(fd, "      end if\n\n");

      var_ptr = var_ptr->next;
   }

   nl = namelists;
   while (nl) {
      if (nl->vtype == INTEGER)
         fortprintf(fd, "      nferr = nf_put_att_int(output_obj %% wr_ncid, NF_GLOBAL, \'%s\', NF_INT, 1, %s)\n", nl->name, nl->name);
      else if (nl->vtype == REAL) {
         fortprintf(fd, "      if (RKIND == 8) then\n", nl->name);
         fortprintf(fd, "         nferr = nf_put_att_double(output_obj %% wr_ncid, NF_GLOBAL, \'%s\', NF_DOUBLE, 1, %s)\n", nl->name, nl->name);
         fortprintf(fd, "      else if (RKIND == 4) then\n", nl->name);
         fortprintf(fd, "         nferr = nf_put_att_real(output_obj %% wr_ncid, NF_GLOBAL, \'%s\', NF_FLOAT, 1, %s)\n", nl->name, nl->name);
         fortprintf(fd, "      end if\n");
      }
      else if (nl->vtype == CHARACTER)
         fortprintf(fd, "      nferr = nf_put_att_text(output_obj %% wr_ncid, NF_GLOBAL, \'%s\', len_trim(%s), trim(%s))\n", nl->name, nl->name, nl->name);
      else if (nl->vtype == LOGICAL) {
         fortprintf(fd, "      if (%s) then\n", nl->name);
         fortprintf(fd, "         nferr = nf_put_att_text(output_obj %% wr_ncid, NF_GLOBAL, \'%s\', 1, \'T\')\n", nl->name);
         fortprintf(fd, "      else\n");
         fortprintf(fd, "         nferr = nf_put_att_text(output_obj %% wr_ncid, NF_GLOBAL, \'%s\', 1, \'F\')\n", nl->name);
         fortprintf(fd, "      end if\n");
      }
      nl = nl->next;
   }

   fclose(fd);   
#endif


   /*
    *  MGD NEW CODE
    */
   fd = fopen("add_output_fields.inc", "w");

   group_ptr = groups;
   while (group_ptr) {
      var_list_ptr = group_ptr->vlist;
      while (var_list_ptr) {
         var_ptr = var_list_ptr->var;

         if (group_ptr->vlist->var->ntime_levs > 1)
            snprintf(struct_deref, 1024, "domain %% blocklist %% %s %% time_levs(1) %% %s", group_ptr->name, group_ptr->name);
         else
            snprintf(struct_deref, 1024, "domain %% blocklist %% %s", group_ptr->name);
         
         if (strncmp(var_ptr->super_array, "-", 1024) != 0) {
            fortprintf(fd, "      if ((%s %% %s %% ioinfo %% output .and. output_obj %% stream == OUTPUT) .or. &\n", struct_deref, var_ptr->super_array);
            fortprintf(fd, "          (%s %% %s %% ioinfo %% restart .and. output_obj %% stream == RESTART) .or. &\n", struct_deref, var_ptr->super_array);
            fortprintf(fd, "          (%s %% %s %% ioinfo %% sfc .and. output_obj %% stream == SFC)) then\n", struct_deref, var_ptr->super_array);
            memcpy(super_array, var_ptr->super_array, 1024);
            fortprintf(fd, "         call MPAS_streamAddField(output_obj %% io_stream, %s %% %s, ierr)\n", struct_deref, super_array);
            while (var_list_ptr && strncmp(super_array, var_list_ptr->var->super_array, 1024) == 0) {
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


   /*
    *  MGD NEW CODE
    */
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

   
#ifdef LEGACY_CODE
   /*
    *  Generate collect and write code
    */
   fd = fopen("io_output_fields.inc", "w");

   group_ptr = groups;
   while (group_ptr) {
      var_list_ptr = group_ptr->vlist;
      while (var_list_ptr) {
         var_ptr = var_list_ptr->var;

         if (group_ptr->vlist->var->ntime_levs > 1)
            snprintf(struct_deref, 1024, "domain %% blocklist %% %s %% time_levs(1) %% %s", group_ptr->name, group_ptr->name);
         else
            snprintf(struct_deref, 1024, "domain %% blocklist %% %s", group_ptr->name);
         
         i = 1;
         dimlist_ptr = var_ptr->dimlist;
         if (var_ptr->vtype == INTEGER) sprintf(vtype, "int"); 
         else if (var_ptr->vtype == REAL) sprintf(vtype, "real"); 
         else if (var_ptr->vtype == CHARACTER) sprintf(vtype, "char"); 
   
         if (strncmp(var_ptr->super_array, "-", 1024) != 0) {
            fortprintf(fd, "      if ((%s %% %s %% ioinfo %% output .and. output_obj %% stream == OUTPUT) .or. &\n", struct_deref, var_ptr->super_array);
            fortprintf(fd, "          (%s %% %s %% ioinfo %% restart .and. output_obj %% stream == RESTART) .or. &\n", struct_deref, var_ptr->super_array);
            fortprintf(fd, "          (%s %% %s %% ioinfo %% sfc .and. output_obj %% stream == SFC)) then\n", struct_deref, var_ptr->super_array);
         }
         else {
            fortprintf(fd, "      if ((%s %% %s %% ioinfo %% output .and. output_obj %% stream == OUTPUT) .or. &\n", struct_deref, var_ptr->name_in_code);
            fortprintf(fd, "          (%s %% %s %% ioinfo %% restart .and. output_obj %% stream == RESTART) .or. &\n", struct_deref, var_ptr->name_in_code);
            fortprintf(fd, "          (%s %% %s %% ioinfo %% sfc .and. output_obj %% stream == SFC)) then\n", struct_deref, var_ptr->name_in_code);
         }
   
         if (var_ptr->ndims > 0) {
            while (dimlist_ptr) {
                  if (i < var_ptr->ndims) {
                     fortprintf(fd, "      %s%id %% ioinfo %% start(%i) = 1\n", vtype, var_ptr->ndims, i);
                     if (dimlist_ptr->dim->constant_value < 0)
                        if (!dimlist_ptr->dim->namelist_defined) fortprintf(fd, "      %s%id %% ioinfo %% count(%i) = domain %% blocklist %% mesh %% %s\n", vtype, var_ptr->ndims, i, dimlist_ptr->dim->name_in_code);
                        else fortprintf(fd, "      %s%id %% ioinfo %% count(%i) = domain %% blocklist %% mesh %% %s\n", vtype, var_ptr->ndims, i, dimlist_ptr->dim->name_in_file);
                     else
                        fortprintf(fd, "      %s%id %% ioinfo %% count(%i) = %s\n", vtype, var_ptr->ndims, i, dimlist_ptr->dim->name_in_code);
                  }
                  else {
                     fortprintf(fd, "      %s%id %% ioinfo %% start(%i) = 1\n", vtype, var_ptr->ndims, i);
                     if (is_derived_dim(dimlist_ptr->dim->name_in_code)) {
                        split_derived_dim_string(dimlist_ptr->dim->name_in_code, &cp1, &cp2);
                        fortprintf(fd, "      %s%id %% ioinfo %% count(%i) = n%sGlobal%s\n", vtype, var_ptr->ndims, i, cp1, cp2);
                        free(cp1);
                        free(cp2);
                     }
                     else if (dimlist_ptr->dim->constant_value > 0) {
                        fortprintf(fd, "      %s%id %% ioinfo %% count(%i) = %s\n", vtype, var_ptr->ndims, i, dimlist_ptr->dim->name_in_code);
                     }
                     else
                        if (!dimlist_ptr->dim->namelist_defined) fortprintf(fd, "      %s%id %% ioinfo %% count(%i) = %sGlobal\n", vtype, var_ptr->ndims, i, dimlist_ptr->dim->name_in_code);
                        else fortprintf(fd, "      %s%id %% ioinfo %% count(%i) = %sGlobal\n", vtype, var_ptr->ndims, i, dimlist_ptr->dim->name_in_file);
                  }
               dimlist_ptr = dimlist_ptr->next;
               i++;
            }
      
            fortprintf(fd, "      allocate(%s%id %% array(", vtype, var_ptr->ndims);
            i = 1;
            dimlist_ptr = var_ptr->dimlist;
      
            if (i < var_ptr->ndims)
               if (dimlist_ptr->dim->constant_value < 0)
                  if (!dimlist_ptr->dim->namelist_defined) fortprintf(fd, "domain %% blocklist %% mesh %% %s", dimlist_ptr->dim->name_in_code);
                  else fortprintf(fd, "domain %% blocklist %% mesh %% %s", dimlist_ptr->dim->name_in_file);
               else
                  fortprintf(fd, "%s", dimlist_ptr->dim->name_in_code);
            else {
               if (is_derived_dim(dimlist_ptr->dim->name_in_code)) {
                  split_derived_dim_string(dimlist_ptr->dim->name_in_code, &cp1, &cp2);
                  fortprintf(fd, "n%sGlobal%s", cp1, cp2);
                  free(cp1);
                  free(cp2);
               }
               else if (dimlist_ptr->dim->constant_value > 0) {
                  fortprintf(fd, "%s", dimlist_ptr->dim->name_in_code);
               }
               else
                  if (!dimlist_ptr->dim->namelist_defined) fortprintf(fd, "%sGlobal", dimlist_ptr->dim->name_in_code);
                  else fortprintf(fd, "%sGlobal", dimlist_ptr->dim->name_in_file);
               lastdim = dimlist_ptr;
            }
            dimlist_ptr = dimlist_ptr->next;
            i++;
            while (dimlist_ptr) {
               if (i < var_ptr->ndims)
                  if (dimlist_ptr->dim->constant_value < 0)
                     if (!dimlist_ptr->dim->namelist_defined) fortprintf(fd, ", domain %% blocklist %% mesh %% %s", dimlist_ptr->dim->name_in_code);
                     else fortprintf(fd, ", domain %% blocklist %% mesh %% %s", dimlist_ptr->dim->name_in_file);
                  else
                     fortprintf(fd, ", %s", dimlist_ptr->dim->name_in_code);
               else {
                  if (is_derived_dim(dimlist_ptr->dim->name_in_code)) {
                     split_derived_dim_string(dimlist_ptr->dim->name_in_code, &cp1, &cp2);
                     fortprintf(fd, ", n%sGlobal%s", cp1, cp2);
                     free(cp1);
                     free(cp2);
                  }
                  else if (dimlist_ptr->dim->constant_value > 0) {
                     fortprintf(fd, ", %s", dimlist_ptr->dim->name_in_code);
                  }
                  else
                     if (!dimlist_ptr->dim->namelist_defined) fortprintf(fd, ", %sGlobal", dimlist_ptr->dim->name_in_code);
                     else fortprintf(fd, ", %sGlobal", dimlist_ptr->dim->name_in_file);
                  lastdim = dimlist_ptr;
               }
               dimlist_ptr = dimlist_ptr->next;
               i++;
            }
            fortprintf(fd, "))\n\n");
   
            if (strncmp(var_ptr->super_array, "-", 1024) != 0) {
               if (var_ptr->ndims > 0) {
                  fortprintf(fd, "      allocate(super_%s%id(", vtype, var_ptr->ndims);
                  i = 1;
                  dimlist_ptr = var_ptr->dimlist;
                  while (dimlist_ptr) {
                     if (dimlist_ptr->dim->constant_value < 0)
                        if (!dimlist_ptr->dim->namelist_defined) fortprintf(fd, "domain %% blocklist %% mesh %% %s", dimlist_ptr->dim->name_in_code);
                        else fortprintf(fd, "domain %% blocklist %% mesh %% %s", dimlist_ptr->dim->name_in_file);
                     else
                        fortprintf(fd, "%s", dimlist_ptr->dim->name_in_code);
      
                     if (i < var_ptr->ndims) fortprintf(fd, ", ");
         
                     dimlist_ptr = dimlist_ptr->next;
                     i++;
                  }
                  fortprintf(fd, "))\n\n");
               }
   
               /* Copy from field to super_ array */
               i = 1;
               dimlist_ptr = var_ptr->dimlist;
               while (i <= var_ptr->ndims) {
                  if (dimlist_ptr->dim->constant_value < 0)
                     if (!dimlist_ptr->dim->namelist_defined) fortprintf(fd, "      do i%i=1,domain %% blocklist %% mesh %% %s\n", i, dimlist_ptr->dim->name_in_code);
                     else fortprintf(fd, "      do i%i=1,domain %% blocklist %% mesh %% %s\n", i, dimlist_ptr->dim->name_in_file);
                  else
                     fortprintf(fd, "      do i%i=1,%s\n", i, dimlist_ptr->dim->name_in_code);
   
                  i++;
                  dimlist_ptr = dimlist_ptr->next;
               }
   
               fortprintf(fd, "         super_%s%id(", vtype, var_ptr->ndims);
               for(i=1; i<=var_ptr->ndims; i++) {
                  fortprintf(fd, "i%i",i);
                  if (i < var_ptr->ndims) fortprintf(fd, ",");
               }
               fortprintf(fd, ") = %s %% %s %% array(", struct_deref, var_ptr->super_array);
               fortprintf(fd, "%s %% index_%s", struct_deref, var_ptr->name_in_code);
               for(i=1; i<=var_ptr->ndims; i++) {
                  fortprintf(fd, ",i%i",i);
               }
               fortprintf(fd, ")\n");
   
               i = 1;
               while (i <= var_ptr->ndims) {
                  fortprintf(fd, "      end do\n");
                  i++;
               }
            }
   
            fortprintf(fd, "      %s%id %% ioinfo %% fieldName = \'%s\'\n", vtype, var_ptr->ndims, var_ptr->name_in_file);
            fortprintf(fd, "      call mpas_dmpar_alltoall_field(domain %% dminfo, &\n");
            if (strncmp(var_ptr->super_array, "-", 1024) != 0)
               fortprintf(fd, "                                super_%s%id, %s%id %% array, &\n", vtype, var_ptr->ndims, vtype, var_ptr->ndims);
            else
               fortprintf(fd, "                                %s %% %s %% array, %s%id %% array, &\n", struct_deref, var_ptr->name_in_code, vtype, var_ptr->ndims);
      
            i = 1;
            dimlist_ptr = var_ptr->dimlist;
            
            if (dimlist_ptr->dim->constant_value < 0)
               if (!dimlist_ptr->dim->namelist_defined) fortprintf(fd, "                                domain %% blocklist %% mesh %% %s", dimlist_ptr->dim->name_in_code);
               else fortprintf(fd, "                                domain %% blocklist %% mesh %% %s", dimlist_ptr->dim->name_in_file);
            else
               fortprintf(fd, "                                %s", dimlist_ptr->dim->name_in_code);
       
            dimlist_ptr = dimlist_ptr->next;
            i++;
            while (dimlist_ptr) {
               if (dimlist_ptr->dim->constant_value < 0)
                  if (!dimlist_ptr->dim->namelist_defined) fortprintf(fd, ", domain %% blocklist %% mesh %% %s", dimlist_ptr->dim->name_in_code);
                  else fortprintf(fd, ", domain %% blocklist %% mesh %% %s", dimlist_ptr->dim->name_in_file);
               else
                  fortprintf(fd, ", %s", dimlist_ptr->dim->name_in_code);
      
               dimlist_ptr = dimlist_ptr->next;
               i++;
            }     
      
            /* 
             *  Need to avoid output_obj in case this is a non-decomposed dimension, in which case 
             *   the send/recv lists are local variables 
             */
            if (strncmp(lastdim->dim->name_in_file,"nCells",1024) != 0 &&
                strncmp(lastdim->dim->name_in_file,"nEdges",1024) != 0 &&
                strncmp(lastdim->dim->name_in_file,"nVertices",1024) != 0 &&
                strncmp(lastdim->dim->name_in_file,"nVertLevels",1024) != 0
               ) {
               if (is_derived_dim(lastdim->dim->name_in_code)) {
                  split_derived_dim_string(lastdim->dim->name_in_code, &cp1, &cp2);
                  fortprintf(fd, ", n%sGlobal%s, &\n", cp1, cp2);
                  fortprintf(fd, "                                send%sList, recv%sList)\n", lastdim->dim->name_in_file+1, lastdim->dim->name_in_file+1);
                  free(cp1);
                  free(cp2);
               }
               else if (lastdim->dim->constant_value > 0) {
                  fortprintf(fd, ", %s, &\n", lastdim->dim->name_in_code);
                  fortprintf(fd, "                                send%sList, recv%sList)\n", lastdim->dim->name_in_file, lastdim->dim->name_in_file);
               }
               else {
                  if (!lastdim->dim->namelist_defined) {
                     fortprintf(fd, ", %sGlobal, &\n", lastdim->dim->name_in_code);
                     fortprintf(fd, "                                send%sList, recv%sList)\n", lastdim->dim->name_in_code+1, lastdim->dim->name_in_code+1);
                  }
                  else {
                     fortprintf(fd, ", %sGlobal, &\n", lastdim->dim->name_in_file);
                     fortprintf(fd, "                                send%sList, recv%sList)\n", lastdim->dim->name_in_file+1, lastdim->dim->name_in_file+1);
                  }
               }
            }
            else {
               if (is_derived_dim(lastdim->dim->name_in_code)) {
                  split_derived_dim_string(lastdim->dim->name_in_code, &cp1, &cp2);
                  fortprintf(fd, ", n%sGlobal%s, &\n", cp1, cp2);
                  fortprintf(fd, "                                output_obj %% send%sList, output_obj %% recv%sList)\n", lastdim->dim->name_in_file+1, lastdim->dim->name_in_file+1);
                  free(cp1);
                  free(cp2);
               }
               else {
                  if (!lastdim->dim->namelist_defined) {
                     fortprintf(fd, ", %sGlobal, &\n", lastdim->dim->name_in_code);
                     fortprintf(fd, "                                output_obj %% send%sList, output_obj %% recv%sList)\n", lastdim->dim->name_in_code+1, lastdim->dim->name_in_code+1);
                  }
                  else {
                     fortprintf(fd, ", %sGlobal, &\n", lastdim->dim->name_in_file);
                     fortprintf(fd, "                                output_obj %% send%sList, output_obj %% recv%sList)\n", lastdim->dim->name_in_file+1, lastdim->dim->name_in_file+1);
                  }
               }
            }
         }
         else {
            fortprintf(fd, "      %s%id %% ioinfo %% fieldName = \'%s\'\n", vtype, var_ptr->ndims, var_ptr->name_in_file);
            fortprintf(fd, "      %s%id %% scalar = %s %% %s %% scalar\n", vtype, var_ptr->ndims, struct_deref, var_ptr->name_in_code);
         }
   
         if (var_ptr->timedim)
            fortprintf(fd, "      if (domain %% dminfo %% my_proc_id == IO_NODE) call mpas_io_output_field_time(output_obj, %s%id)\n", vtype, var_ptr->ndims);
         else
            fortprintf(fd, "      if (domain %% dminfo %% my_proc_id == IO_NODE) call mpas_io_output_field(output_obj, %s%id)\n", vtype, var_ptr->ndims);
         if (var_ptr->ndims > 0) {
            fortprintf(fd, "      deallocate(%s%id %% array)\n", vtype, var_ptr->ndims);
            if (strncmp(var_ptr->super_array, "-", 1024) != 0)
               fortprintf(fd, "      deallocate(super_%s%id)\n", vtype, var_ptr->ndims);
         }
         fortprintf(fd, "      end if\n\n");
   
         var_list_ptr = var_list_ptr->next;
      }
      group_ptr = group_ptr->next;
   }

   fclose(fd);
   
   
   /*
    *  Generate code to write 0d, 1d, 2d, 3d time-invariant fields
    */
   for(j=0; j<3; j++) {
      for(i=0; i<=3; i++) {
         if (j == 0) {
            sprintf(fname, "output_field%idinteger.inc", i);
            ivtype = INTEGER;
         }
         else if (j == 1) {
            sprintf(fname, "output_field%idreal.inc", i);
            ivtype = REAL;
         }
         else if (j == 2) {
            sprintf(fname, "output_field%idchar.inc", i);
            ivtype = CHARACTER;
         }
         fd = fopen(fname, "w");
   
         var_ptr = vars;
         while (var_ptr && (var_ptr->ndims != i || var_ptr->vtype != ivtype || var_ptr->timedim)) var_ptr = var_ptr->next;
         if (var_ptr) {
            fortprintf(fd, "      if (trim(field %% ioinfo %% fieldName) == \'%s\') then\n", var_ptr->name_in_file);
            fortprintf(fd, "         varID = output_obj %% wrVarID%s\n", var_ptr->name_in_file);
            var_ptr = var_ptr->next;
            while (var_ptr) {
               if (var_ptr->ndims == i && var_ptr->vtype == ivtype && !var_ptr->timedim) {
                  fortprintf(fd, "      else if (trim(field %% ioinfo %% fieldName) == \'%s\') then\n", var_ptr->name_in_file);
                  fortprintf(fd, "         varID = output_obj %% wrVarID%s\n", var_ptr->name_in_file);
               }
               var_ptr = var_ptr->next;
            }
            fortprintf(fd, "      end if\n");
         }
      
         fclose(fd);
      } 
   } 

   
   /*
    *  Generate code to write 0d, 1d, 2d, 3d time-varying fields
    */
   for(j=0; j<2; j++) {
      for(i=0; i<=3; i++) {
         if (j == 0) {
            sprintf(fname, "output_field%idinteger_time.inc", i);
            ivtype = INTEGER;
         }
         else {
            sprintf(fname, "output_field%idreal_time.inc", i);
            ivtype = REAL;
         }
         fd = fopen(fname, "w");
   
         var_ptr = vars;
         while (var_ptr && (var_ptr->ndims != i || var_ptr->vtype != ivtype || !var_ptr->timedim)) var_ptr = var_ptr->next;
         if (var_ptr) {
            fortprintf(fd, "      if (trim(field %% ioinfo %% fieldName) == \'%s\') then\n", var_ptr->name_in_file);
            fortprintf(fd, "         varID = output_obj %% wrVarID%s\n", var_ptr->name_in_file);
            var_ptr = var_ptr->next;
            while (var_ptr) {
               if (var_ptr->ndims == i && var_ptr->vtype == ivtype && var_ptr->timedim) {
                  fortprintf(fd, "      else if (trim(field %% ioinfo %% fieldName) == \'%s\') then\n", var_ptr->name_in_file);
                  fortprintf(fd, "         varID = output_obj %% wrVarID%s\n", var_ptr->name_in_file);
               }
               var_ptr = var_ptr->next;
            }
            fortprintf(fd, "      end if\n");
         }
      
	 fclose(fd);
      }
   }

   
   /*
    *  Generate code to write 0d and 1d character time-varying fields
    */
   for(i=0; i<=1; i++) {
      sprintf(fname, "output_field%idchar_time.inc", i);
      fd = fopen(fname, "w");

      var_ptr = vars;
      while (var_ptr && (var_ptr->ndims != i || var_ptr->vtype != CHARACTER || !var_ptr->timedim)) var_ptr = var_ptr->next;
      if (var_ptr) {
         fortprintf(fd, "      if (trim(field %% ioinfo %% fieldName) == \'%s\') then\n", var_ptr->name_in_file);
         fortprintf(fd, "         varID = output_obj %% wrVarID%s\n", var_ptr->name_in_file);
         var_ptr = var_ptr->next;
         while (var_ptr) {
            if (var_ptr->ndims == i && var_ptr->vtype == CHARACTER && var_ptr->timedim) {
               fortprintf(fd, "      else if (trim(field %% ioinfo %% fieldName) == \'%s\') then\n", var_ptr->name_in_file);
               fortprintf(fd, "         varID = output_obj %% wrVarID%s\n", var_ptr->name_in_file);
            }
            var_ptr = var_ptr->next;
         }
         fortprintf(fd, "      end if\n");
      }
   
      fclose(fd);
   }
#endif
   
}
