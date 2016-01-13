// The DyCoreToGlimmer class provides methods to move Glimmer Fortran data to C++ structures
// for access by the Chombo-based FELIX model.  The structure names and structure member
// names mostly correspond to derived types defined in Glimmer.  In general, pointers to
// the Fortran data arrays are used, rather than copies of these arrays.  This saves space
// and reduces the steps needed to update the Glimmer data between calls to the FELIX
// ice sheet modeling program.  Methods are provided to set these array pointers, and 
// copy array dimension information.  Objects of this class are accessed by extern C
// routines in felix_to_glim_extern.cpp, and by the FELIX front end.  DMR--5/24/10

#include "FelixToGlimmer.H"


using namespace std;


int
FelixToGlimmer::initDyCore(const char * input_fname)
{

  //  long * dimInfo;

  cout << "In FELIX initDyCore" << endl;
  //  dimInfo = this -> getLongVar("dimInfo","geometry");
    
       
  //  cout << "DimInfo in initDyCore: " << endl;
  //  for (i=0;i<10;i++) cout << dimInfo[i] << " ";     
  //  cout << "In FELIX initDyCore, calling felix_driver_inin:" << endl;
  felix_driver_init(2,0,this,input_fname);
  return 0; // ought to make sensible use of this.

}

// updates cur_time_yr to match time update in dycore
int
FelixToGlimmer::runDyCore(float& cur_time_yr, const float time_inc_yr)
{
  cout << "In FELIX runDyCore" << endl;
  felix_driver_run(this,cur_time_yr,time_inc_yr);
  return 0; // ought to make sensible use of this.
}

int
FelixToGlimmer::deleteDyCore()
{
  felix_driver_finalize(this -> getDyCoreIndex());
  return  0; // ought to make sensible use of this.
}
  
//int storeFelixObject(AmrIce bisicles_object)
//{}

//AmrIce retrieveFelixObject()
//{}
