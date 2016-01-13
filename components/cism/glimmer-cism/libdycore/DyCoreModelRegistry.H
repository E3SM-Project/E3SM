// DyCoreModelRegistry is used to store multiple dynamic core models
// without using global variables, and allowing retrieval by an
// integer model index.  --DMR 5/24/10

#include "DyCoreToGlimmer.H"
//IK, 8/6/13: added some ifdefs here to allow for other external dycores than just bisicles

#ifdef CISM_HAS_BISICLES
#include "BISICLES/BisiclesToGlimmer.H"
#endif
//#include "Ymir/YmirToGlimmer.H"
//IK, 8/6/13: added the following for FELIX as external dycore
#ifdef CISM_HAS_FELIX
#include "FELIX/FelixToGlimmer.H"
#endif 

#ifndef DYCORE_MODEL_REGISTRY
#define DYCORE_MODEL_REGISTRY

#define DYCORE_MODEL_COUNT 10

class DyCoreModelRegistry
{
  private:
  
    struct RegistryEntry {
      DyCoreToGlimmer *  dycore_to_glimmer;
      int dycore_type; // 0=BISICLES, 1=Ymir, 2=FELIX
      int my_reg_index;
      int dycore_present;
    } entry[DYCORE_MODEL_COUNT];

  int cur_model_count;


  public: 
  
  DyCoreModelRegistry();
  int ClearRegistryEntries();
  int ClearRegistryEntry(int index);
  DyCoreToGlimmer * getDyCoreToGlimmerByIndex(int index);
  //int * getDyCoreToGlimmerByIndex(int index);
  int getModelCount();
  int incModelCount();
  int setDyCoreByType(int index,int dyncore_type);
  int setRegistryIndex(int index);
  int getRegistryIndex(int index); 
};

#endif
