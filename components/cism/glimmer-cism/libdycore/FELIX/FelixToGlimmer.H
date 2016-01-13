// The DyCoreToGlimmer class provides methods to move Glimmer Fortran data to C++ structures
// for access by the Chombo-based FELIX model.  The structure names and structure member
// names mostly correspond to derived types defined in Glimmer.  In general, pointers to
// the Fortran data arrays are used, rather than copies of these arrays.  This saves space
// and reduces the steps needed to update the Glimmer data between calls to the FELIX
// ice sheet modeling program.  Methods are provided to set these array pointers, and 
// copy array dimension information.  Objects of this class are accessed by extern C
// routines in felix_to_glim_extern.cpp, and by the FELIX front end.  DMR--5/24/10

#ifndef FELIXTOGLIMMER
#define FELIXTOGLIMMER


#include <stdio.h>
#include <iostream>
#include "/home/ikalash/Desktop/clean/Albany2/src/FELIX/interface_with_cism/felix_driver.H"
#include "../DyCoreToGlimmer.H"


class FelixToGlimmer : public DyCoreToGlimmer
{
  private:
  // AmrIce bisicles_object;

  public: 
  
  //DynCoreToGlimmer FelixToGlimmer();
  int initDyCore(const char * input_fname);
  // cur_time_yr is updated in place as solution is evolved
  int runDyCore(float& cur_time_yr, const float time_inc_yr);
  int deleteDyCore();

};

#endif
