// DyCoreModelRegistry is used to store multiple dynamic core models
// without using global variables, and allowing retrieval by an
// integer model index.  --DMR 5/24/10

//#include "DyCoreToGlimmer.H"
#include "BISICLES/BisiclesToGlimmer.H"
#include "Ymir/YmirToGlimmer.H"

#ifndef DYCORE_MODEL_REGISTRY
#define DYCORE_MODEL_REGISTRY

#define DYCORE_MODEL_COUNT 10

class DyCoreModelRegistry
{
  private:
  
    struct RegistryEntry {
      DyCoreToGlimmer *  dycore_to_glimmer;
      int dycore_type; // 0=BISICLES, 1=Ymir
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
