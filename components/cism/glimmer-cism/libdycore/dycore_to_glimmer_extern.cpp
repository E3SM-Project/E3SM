// dyncore_to_glimmer_extern.cpp contains the extern C routines used to provide an interface
// between the Fortran based Glimmer and the C++ based BISICLES.  These routines access the
// dycore_registry routine to create and use DyCoreToGlimmer class objects.  Since the registry
// can contain multiple DyCoreToGlimmer objects, it allows multiple DYCORE models to be
// instantiated and used by Glimmer. DMR--5/24/10

#include <malloc.h>
#include <cstring>
#include <stdio.h>
#include <iostream>
#include "DyCoreModelRegistry.H"


extern "C" {
  void dycore_init_registry_();
  void dycore_reset_registry_();
  void dycore_get_new_model_(int * dycore_type,int * index,int * error_code);
  void dycore_init_model_(int * dycore_type,int * index,char * input_fname,int * error_code);
  void dycore_run_model_(int * model_index, float * cur_time_yr, float * time_inc_yr);
  void dycore_delete_model_(int * dycore_model_index);

  void dycore_set_ptr_double_var_(double *var, char *var_name_in,
			          char *struct_name_in, int *model_index);
  void dycore_set_ptr_long_var_(long *var, char *var_name,
			        char *struct_name, int *model_index);

  void dycore_copy_in_double_var_(double *var, char *var_name, char *struct_name,
				  long *dim_info, int *model_index);
  void dycore_copy_in_long_var_(long *var, char *var_name, char *struct_name,
                                long *dim_info, int *model_index);  

  void dycore_test_vel_input_(int *model_index,double *vel,char *var_name,int *dim_info);
  void dycore_test_vel_output_(int *model_index,double *vel,char *var_name,int *dim_info);

}

using namespace std;


// dycore_registry uses the DyCoreModelRegistry class to create a registry of DyCoreToGlimmer
// objects that are interface instances between DYCORE and Glimmer.  This is the only
// routine in this file that accesses a DyCoreModelRegistry object.  It is not accessed
// directly from  Fortran, but used by the extern routines in this file.
int dycore_registry(int init,int get_model_by_index,int * model_index, 
	   DyCoreToGlimmer ** dycore_to_glimmer_ptr,int dycore_type,int clear_entry)
{
  // this declaration initializes the registry, when dycore_registry
  // is first called:
  static DyCoreModelRegistry dmr;

  if (init == 1){
    cout << "Initializing Dycore Model Registry" << endl;
    return(0);
  }
  if (init == -1){
    dmr.ClearRegistryEntries();
    cout << "Cleared Dycore Model Registry" << endl;
    return(0);
  }
  if (clear_entry > 0) {
    cout << "Calling ClearRegistryEntry, for entry:  " << clear_entry << endl;
    dmr.ClearRegistryEntry(clear_entry);
    return(0);
  } 

  if (get_model_by_index == 1) {
    if (*model_index == -1) {
      // if model_index=-1, initialize a new registry entry and
      // obtain a new model index:
      *model_index = dmr.getModelCount() + 1;

      // init a dycore interface object, and add it to the registry:
      dmr.setDyCoreByType(*model_index,dycore_type);

      dmr.setRegistryIndex(*model_index);
     
      if (*model_index > DYCORE_MODEL_COUNT) {
        cout << "Error, exceeded DYCORE Registry limit of " << 
	DYCORE_MODEL_COUNT << endl;
        return(-1);
      }
      dmr.incModelCount();
      return(0);
    }
    // get pointer to DyCoreToGlimmer object from registry:    
    *dycore_to_glimmer_ptr = dmr.getDyCoreToGlimmerByIndex(*model_index);

    return(dmr.getRegistryIndex(*model_index));
  }
  return(0);
}

void dycore_init_registry_()
{
  DyCoreToGlimmer * dummy_dtg;
  int init_registry = 1;

  // initialize a registry if dycore model interfaces:
  dycore_registry(init_registry,0,0,&dummy_dtg,0,0);
}

void dycore_reset_registry_()
{
  DyCoreToGlimmer * dummy_dtg;
  int init_registry = -1;  // set init_registry to clear registry

  // initialize a registry if dycore model interfaces:
  dycore_registry(init_registry,0,0,&dummy_dtg,0,0);
}


void dycore_get_new_model_(int * dycore_type,int * index,int * error_code)
{
  DyCoreToGlimmer * dtg;
  int model_index=-1;

  // cout << "In dycore_get_new_model_ , dycore_type = " << *dycore_type << endl;

  // use *model_index=-1 to initialize a new registry entry:
  *error_code = dycore_registry(0,1,&model_index,&dtg,*dycore_type,0); 
  *index = model_index; 
}

void dycore_init_model_(int * dycore_type,int * model_index,char * input_fname,int * error_code)
{
  DyCoreToGlimmer * dtg;

  // cout << "In dycore_init_model_ , dycore_type = " << *dycore_type << endl;

  dycore_registry(0,1,model_index,&dtg,-1,0);
  
  dtg -> setDyCoreType(*dycore_type);
  dtg -> setDyCoreIndex(*model_index); 
  dtg -> initDyCore(input_fname);
}

void dycore_run_model_(int * model_index, float * cur_time_yr, float * time_inc_yr)
{
  DyCoreToGlimmer * dtg;
  
  dycore_registry(0,1,model_index,&dtg,-1,0); 

  //cout << "In dycore_run_model, model_index = " << *model_index << endl;
  //cout << "In drm, cur_time, time_inc = " << *cur_time_yr << "  " << *time_inc_yr << endl;

  dtg -> runDyCore(*cur_time_yr,*time_inc_yr);
}

void dycore_delete_model_(int * model_index)
{
  DyCoreToGlimmer * dtg;
  int clear_entry;

  clear_entry = *model_index;
  dycore_registry(0,1,model_index,&dtg,-1,clear_entry);
  //  reg_index = dycore_registry(0,1,model_index,&dtg,-1);
  //  dtg -> deleteDyCore();   
}

void dycore_set_ptr_double_var_(double *var, char *var_name,
                           char *struct_name, int *model_index)
{
  DyCoreToGlimmer * dtg;

  dycore_registry(0,1,model_index,&dtg,-1,0);
  dtg -> setDoubleVar(var,var_name,struct_name);
}

void dycore_set_ptr_long_var_(long *var, char *var_name, int *var_name_len,
                           char *struct_name, int *struct_name_len, int *model_index)
{
  DyCoreToGlimmer * dtg;

  // cout << "var_name::" << var_name << "::" << endl;
 
  dycore_registry(0,1,model_index,&dtg,-1,0);
  dtg -> setLongVar(var,var_name,struct_name);
}

void dycore_copy_in_double_var_(double *var, char *var_name, char *struct_name,
                                long *dim_info, int *model_index)
{
  DyCoreToGlimmer * dtg;

  // cout << "In copy_in_double_var, var_name::" << var_name << "::" << endl;
 std::cout  << " dycore_copy_in_double_var_ " << var_name 
	    << " = " << *var << std::endl;
  dycore_registry(0,1,model_index,&dtg,-1,0);
  dtg -> copyInDoubleVar(var,var_name,struct_name,dim_info);
}

void dycore_copy_in_long_var_(long *var, char *var_name, char *struct_name,
                              long *dim_info, int *model_index)
{
  DyCoreToGlimmer * dtg;

  //  cout << "In copy_long_var" << endl;
  //cout << "struct_name::" << struct_name << "::" << endl;
 
  dycore_registry(0,1,model_index,&dtg,-1,0);
  dtg -> copyInLongVar(var,var_name,struct_name,dim_info);
}


void dycore_test_vel_input_(int *model_index,double *vel,char *var_name,
                              int * dim_info)
{
  int i, reg_index;
  //  double test_array[10];
  DyCoreToGlimmer  * dtg;
  //  double * var;

  // cout << "test_vel_in, Calling dycore_registry" << endl;

  reg_index = dycore_registry(0,1,model_index,&dtg,-1,0);

  cout << "test_vel_in model_index compare: " << *model_index << "  " 
                                              << reg_index << endl; 

  //  if (*model_index == 1) (*dtg).set_velocity_data(vel,"uvel",dim_info);
  //  if (*model_index == 2) (*dtg).set_velocity_data(vel,"vvel",dim_info);

  cout << "In vel input, var = " << var_name << ": ";
  for (i=0;i<14;i++) cout << vel[i] << " ";
  cout << endl;
}

void dycore_test_vel_output_(int *model_index,double *vel,char *var_name,
                               int * dim_info)
{
  int i;
  //  double test_array[10];
  DyCoreToGlimmer  * dtg;
  double * var;

  cout << "In test output, Model Index: " << *model_index << endl;
 
  cout << "output: my_reg_index: " << dycore_registry(0,1,model_index,&dtg,-1,0) << endl;

  if (*model_index == 1) {
    //    (*dtg).get_velocity_data(&var,"uvel",dim_info);
  }
   
  if (*model_index == 2) { 
    //    (*dtg).get_velocity_data(&var,"vvel",dim_info);
  }
  cout << "In vel output, var = " << var_name << ": ";
  for (i=0;i<14;i++) cout << var[i] << " ";
  cout << endl;

}
