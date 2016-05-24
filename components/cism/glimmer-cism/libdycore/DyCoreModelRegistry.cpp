// DyCoreModelRegistry is used to store multiple dynamic core models
// without using global variables, and allowing retrieval by an
// integer model index.  --DMR 5/24/10

// 4/9/12 (DMR) commented out Ymir object creation because of a build
// problem on hopper.nersc.gov

#include <stdio.h>
#include <iostream>
#include "DyCoreModelRegistry.H"

using namespace std;

// RegistryEntry entry[BISICLES_MODEL_COUNT];

DyCoreModelRegistry::DyCoreModelRegistry()
{
  cur_model_count = 0;
}

int
DyCoreModelRegistry::ClearRegistryEntries()
{
  int i;
 
  for (i=0;i<DYCORE_MODEL_COUNT;i++) {
    if (entry[i].dycore_present == 1) {
      ClearRegistryEntry(i);
    }
  }
  cur_model_count = 0;
  return(0);
}

int
DyCoreModelRegistry::ClearRegistryEntry(int index)
{
  DyCoreToGlimmer * dtg;

  if (entry[index].dycore_present == 1) {
    cout << "Clearing Registry Entry: " << index << endl;
    dtg = entry[index].dycore_to_glimmer;
    dtg -> deleteDyCore();
    delete entry[index].dycore_to_glimmer;
  }
  entry[index].dycore_present = 0;
  return(0);
}

DyCoreToGlimmer *
DyCoreModelRegistry::getDyCoreToGlimmerByIndex(int index)
{
  
  // cout << index << " Registry entry dycore type: " << entry[index].dycore_type << endl;
  return((DyCoreToGlimmer *) entry[index].dycore_to_glimmer);
}

int
DyCoreModelRegistry::setDyCoreByType(int index,int dycore_type)
{
  entry[index].dycore_type = dycore_type;
  entry[index].dycore_present = 1;

  switch (entry[index].dycore_type) {
    case 0: 
       entry[index].dycore_to_glimmer = new BisiclesToGlimmer;
       break;
    case 1:
      //entry[index].dycore_to_glimmer = new YmirToGlimmer;
       break;
    default: entry[index].dycore_to_glimmer = new BisiclesToGlimmer;
  }
  return(0);
}


int
DyCoreModelRegistry::getModelCount()
{
  return(cur_model_count);
}

int
DyCoreModelRegistry::incModelCount()
{
  cur_model_count++;
  return(0);
}

int
DyCoreModelRegistry::setRegistryIndex(int index)
{
  entry[index].my_reg_index = index;
  return(0);
}

int
DyCoreModelRegistry::getRegistryIndex(int index)
{
  return(entry[index].my_reg_index);
}
