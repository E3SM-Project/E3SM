// The DyCoreToGlimmer class provides methods to move Glimmer Fortran data to C++ structures
// for access by the Chombo-based BISICLES model.  The structure names and structure member
// names mostly correspond to derived types defined in Glimmer.  In general, pointers to
// the Fortran data arrays are used, rather than copies of these arrays.  This saves space
// and reduces the steps needed to update the Glimmer data between calls to the BISICLES
// ice sheet modeling program.  Methods are provided to set these array pointers, and 
// copy array dimension information.  Objects of this class are accessed by extern C
// routines in dycore_to_glimmer_extern.cpp, and by the BISICLES front end.  DMR--5/24/10

#include "DyCoreToGlimmer.H"
#include <cstring>

using namespace std;

//DyCoreToGlimmer::DyCoreToGlimmer(int dycore_type)
DyCoreToGlimmer::DyCoreToGlimmer()
{
}

DyCoreToGlimmer::~DyCoreToGlimmer()
{
  // cout << "Init DyCoreToGlimmer" << endl;
}

int
DyCoreToGlimmer::setDoubleVar(double *var, const char *var_name,  const char *struct_name)
{
  //cout << "struct_name::" << struct_name << "::" << endl;

  if (strcmp(struct_name,"geometry") == 0) {
    if (strcmp(var_name,"thck") == 0) geometry.thck = var;
    if (strcmp(var_name,"topg") == 0) geometry.topg = var;
    if (strcmp(var_name,"usrf") == 0) geometry.usrf = var;
  }

  if (strcmp(struct_name,"velocity") == 0) {
    if (strcmp(var_name,"uvel") == 0) velocity.uvel = var; 
    if (strcmp(var_name,"vvel") == 0) velocity.vvel = var;
    if (strcmp(var_name,"wvel") == 0) velocity.wvel = var;
    if (strcmp(var_name,"btrc") == 0) velocity.btrc = var;
    //   if (strcmp(var_name,"wgrd") == 0) velocity.wgrd = var;
    //cout << "Set velocity var, " << var_name << endl;
  }

 

  if (strcmp(struct_name,"temper") == 0) {
    if (strcmp(var_name,"temp") == 0) temper.temp = var;
    if (strcmp(var_name,"bheatflx") == 0) temper.bheatflx = var;
    if (strcmp(var_name,"bmlt") == 0) temper.bmlt = var;
  }  
  
  if (strcmp(struct_name,"climate") == 0) {
    if (strcmp(var_name,"acab") == 0) climate.acab = var;
    if (strcmp(var_name,"acab_tavg") == 0) climate.acab_tavg = var;
    if (strcmp(var_name,"calving") == 0) climate.calving = var;
  }

  return(0);
}

double *
DyCoreToGlimmer::getDoubleVar(const char *var_name, const  char *struct_name)
{
  
  double * var=0;

  //cout << "struct_name::" << struct_name << "::" << endl;

  if (strcmp(struct_name,"geometry") == 0) {
    if (strcmp(var_name,"thck") == 0)
    {
      return(geometry.thck);
    }
    if (strcmp(var_name,"topg") == 0)
    {
      return(geometry.topg);
    }
    if (strcmp(var_name,"usrf") == 0)
    {
      return(geometry.usrf);
    }
  }
  if (strcmp(struct_name,"numerics") == 0) {
    if (strcmp(var_name,"dew") == 0) return(numerics.dew);
    if (strcmp(var_name,"dns") == 0) return(numerics.dns);
  }
  if (strcmp(struct_name,"velocity") == 0) {
    if (strcmp(var_name,"btrc") == 0) return (velocity.btrc);
    if (strcmp(var_name,"uvel") == 0) return (velocity.uvel);
    if (strcmp(var_name,"vvel") == 0) return (velocity.vvel);
    if (strcmp(var_name,"wvel") == 0) return (velocity.wvel);

    //cout << "Set velocity var, " << var_name << endl;
  }

  if (strcmp(struct_name,"temper") == 0) {
    if (strcmp(var_name,"temp") == 0) var = temper.temp;
    if (strcmp(var_name,"bheatflx") == 0) var = temper.bheatflx;
    if (strcmp(var_name,"bmlt") == 0) var = temper.bmlt;
  }

  if (strcmp(struct_name,"climate") == 0) {
    if (strcmp(var_name,"acab") == 0) var = climate.acab;
  }

  return(var);
}


int
DyCoreToGlimmer::setLongVar(long * var,  const char *var_name,  const char *struct_name)
{
  // cout << "struct_name::" << struct_name << "::" << endl;

  if (strcmp(struct_name,"geometry") == 0) {
    if (strcmp(var_name,"dimInfo") == 0) geometry.dimInfo = var;
  }
  if (strcmp(struct_name,"velocity") == 0) {
    //cout << "Set velocity var, " << var_name << endl;
  }
  return(0);
}

long * 
DyCoreToGlimmer::getLongVar( const char *var_name,  const char *struct_name)
{
  long * var;

  if (strcmp(struct_name,"geometry") == 0) {
    if (strcmp(var_name,"dimInfo") == 0) var = geometry.dimInfo;
    if (strcmp(var_name,"ewlb") == 0) var = geometry.ewlb;
    if (strcmp(var_name,"ewub") == 0) var = geometry.ewub;
    if (strcmp(var_name,"nslb") == 0) var = geometry.nslb;
    if (strcmp(var_name,"nsub") == 0) var = geometry.nsub;
  }
  if (strcmp(struct_name,"mpi_vars") == 0) {
    if (strcmp(var_name,"communicator") == 0) var = mpi_vars.communicator;
    if (strcmp(var_name,"process_count") == 0) var = mpi_vars.process_count;
    if (strcmp(var_name,"my_rank") == 0) var = mpi_vars.my_rank;
  }

  if (strcmp(struct_name,"velocity") == 0) {
    if (strcmp(var_name,"dimInfo") == 0) var = velocity.dimInfo;
  }

  if (strcmp(struct_name,"climate") == 0) {
    if (strcmp(var_name,"dimInfo") == 0) var = climate.dimInfo;
  }

  return(var);
}


int
DyCoreToGlimmer::copyInDoubleVar( const double *var,  const char *var_name, 
				  const char *struct_name, const long *var_dim_info)
{
  long elem_count=1;
  long i;
  
  // std::cout  << "copyInDoubleVar " << var_name << " = " << *var << std::endl;

  for (i=1;i<=var_dim_info[0];i++) elem_count *= var_dim_info[i];
    
  //cout << "struct_name::" << struct_name << "::" << endl;
  if (strcmp(struct_name,"geometry") == 0) {
    if (strcmp(var_name,"dimInfo") == 0) {
      
    }
  }
 
  if (strcmp(struct_name,"velocity") == 0) {

  }
 
  if (strcmp(struct_name,"numerics") == 0) {
    if (strcmp(var_name,"dew") == 0) {
      numerics.dew = new double[elem_count];
      for (i=0;i<elem_count;i++) numerics.dew[i] = var[i];
      //cout << "Copy dew. dew, elem_count = "  << *numerics.dew << "  " << elem_count << endl;
    }  
    if (strcmp(var_name,"dns") == 0) {
      numerics.dns = new double[elem_count];
      for (i=0;i<elem_count;i++) numerics.dns[i] = var[i];
      //cout << "Copy dns. dns, elem_count = "  << *numerics.dns << "  " << elem_count << endl;
    }
  }
  if (strcmp(struct_name,"climate") == 0) {
    if (strcmp(var_name,"eus") == 0) {
      climate.eus = new double[elem_count];
      for (i=0;i<elem_count;i++) climate.eus[i] = var[i];
    }
  }
  return(0);
}


int
DyCoreToGlimmer::copyInLongVar(const long *var, const char *var_name, 
			       const char *struct_name, const long *var_dim_info)
{
  long elem_count=1;
  long i;

  //cout << "DycoreToGlimmer copy long, size = " << var_dim_info[0] << endl; 
  for (i=1;i<=var_dim_info[0];i++) elem_count *= var_dim_info[i];
    
  //cout << "DycoreToGlimmer -- struct_name::" << struct_name << "::" << sizeof(long) << endl;
  if (strcmp(struct_name,"velocity") == 0) {
    if (strcmp(var_name,"dimInfo") == 0) {
      velocity.dimInfo = new long[elem_count];
      for (i=0;i<elem_count;i++) velocity.dimInfo[i] = var[i];
    }  
  }
  if (strcmp(struct_name,"geometry") == 0) {
    if (strcmp(var_name,"dimInfo") == 0) {
      geometry.dimInfo = new long[elem_count];
      for (i=0;i<elem_count;i++) geometry.dimInfo[i] = var[i];
    }  
    if (strcmp(var_name,"ewlb") == 0) {
      geometry.ewlb = new long[elem_count];
      for (i=0;i<elem_count;i++) geometry.ewlb[i] = var[i];
    }  
    if (strcmp(var_name,"ewub") == 0) {
      geometry.ewub = new long[elem_count];
      for (i=0;i<elem_count;i++) geometry.ewub[i] = var[i];
    }  
    if (strcmp(var_name,"nslb") == 0) {
      geometry.nslb = new long[elem_count];
      for (i=0;i<elem_count;i++) geometry.nslb[i] = var[i];
    }  
    if (strcmp(var_name,"nsub") == 0) {
      geometry.nsub = new long[elem_count];
      for (i=0;i<elem_count;i++) geometry.nsub[i] = var[i];
    }  
  }
  if (strcmp(struct_name,"climate") == 0) {
    if (strcmp(var_name,"dimInfo") == 0) {
      climate.dimInfo = new long[elem_count];
      for (i=0;i<elem_count;i++) climate.dimInfo[i] = var[i];
    }  
  }
  if (strcmp(struct_name,"mpi_vars") == 0) {
    if (strcmp(var_name,"communicator") == 0) {
      mpi_vars.communicator = new long[elem_count];
      for (i=0;i<elem_count;i++) mpi_vars.communicator[i] = var[i];
    }
  }
  if (strcmp(struct_name,"mpi_vars") == 0) {
    if (strcmp(var_name,"process_count") == 0) {
      mpi_vars.process_count = new long[elem_count];
      for (i=0;i<elem_count;i++) mpi_vars.process_count[i] = var[i];
    }
  }
  if (strcmp(struct_name,"mpi_vars") == 0) {
    if (strcmp(var_name,"my_rank") == 0) {
      mpi_vars.my_rank = new long[elem_count];
      for (i=0;i<elem_count;i++) mpi_vars.my_rank[i] = var[i];
    }
  }

  return(0);
}


int 
DyCoreToGlimmer::initDyCore(const char * dycore_fname)
{
  cout << "In DycoreToGlimmer::initDyCore" << endl;
  return(0);
}

int
DyCoreToGlimmer::runDyCore(float& cur_time_yr, const float time_inc_yr)
{
  return(0);
}

int 
DyCoreToGlimmer::deleteDyCore()
{
  return(0);
}

int 
DyCoreToGlimmer::setDyCoreType(const int dycore_type)
{
  //cout << "in set type, type = " << dycore_type << endl;
  dycore_info.dycore_type = dycore_type;
  return(0);
}

int
DyCoreToGlimmer::getDyCoreType()
{
  return(dycore_info.dycore_type);
}

int 
DyCoreToGlimmer::setDyCoreIndex(const int dycore_index)
{
  //cout << "in set index, index = " << dycore_index << endl;
  dycore_info.dycore_index = dycore_index;
  return(0);
}

int
DyCoreToGlimmer::getDyCoreIndex()
{
  return(dycore_info.dycore_index);
}
