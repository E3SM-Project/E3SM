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
  // initialize all pointerrs to null
  // geometry...
  geometry.thck = NULL;
  geometry.usrf = NULL;
  geometry.lsrf = NULL;
  geometry.topg = NULL;
  geometry.floating_mask = NULL;
  geometry.ice_mask = NULL;
  geometry.lower_cell_loc = NULL;
  geometry.lower_cell_temp = NULL;
  geometry.dimInfo = NULL;
  geometry.ewlb = NULL;
  geometry.ewub = NULL;
  geometry.nslb = NULL;
  geometry.nsub = NULL;
  geometry.nhalo = NULL;
  
  // velocity
  velocity.uvel = NULL; //output
  velocity.vvel = NULL; //output
  velocity.wvel = NULL;
  velocity.wgrd = NULL;
  velocity.btrc = NULL; // basal traction coefficient
  velocity.dimInfo = NULL;

  // temper
  temper.temp = NULL; // Three-dimensional temperature field.
  temper.bheatflx = NULL; // basal heat flux (2D)
  temper.bmlt = NULL;  // Basal melt-rate

  // numerics
  numerics.tstart = NULL; // start time
  numerics.tend = NULL; // end time
  numerics.time = NULL; // current time

  numerics.dew = NULL; // ew cell size
  numerics.dns = NULL; // ns cell size

  // constants are doubles, not pointers. Set to nonphysical values
  double bogusVal = -1.2345e10;
  constants.seconds_per_year = bogusVal;
  constants.gravity = bogusVal;
  constants.rho_ice = bogusVal; 
  constants.rho_seawater = bogusVal;  
  constants.therm_diffusivity_ice = bogusVal;
  constants.heat_capacity_ice = bogusVal;


  // climate
  climate.acab = NULL; // Annual mass balance.
  climate.acab_tavg = NULL; // Annual mass balance (time average)
  climate.calving = NULL; // Calving flux (scaled as mass balance, thickness,)
  climate.dimInfo = NULL;
  climate.eus = NULL; // eustatic sea level

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
    else if (strcmp(var_name,"topg") == 0) geometry.topg = var;
    else if (strcmp(var_name,"usrf") == 0) geometry.usrf = var;
    else if (strcmp(var_name,"lsrf") == 0) geometry.lsrf = var;
    else if (strcmp(var_name,"floating_mask") == 0) geometry.floating_mask = var;
    else if (strcmp(var_name,"ice_mask") == 0) geometry.ice_mask = var;
    else if (strcmp(var_name,"lower_cell_loc") == 0) geometry.lower_cell_loc = var;
    else if (strcmp(var_name,"lower_cell_temp") == 0) geometry.lower_cell_temp = var;
    else cerr << "unknown variable type = " << struct_name
              << "." << var_name << " undefined!" << endl;
  }

  else if (strcmp(struct_name,"velocity") == 0) {
    if (strcmp(var_name,"uvel") == 0) velocity.uvel = var; 
    else if (strcmp(var_name,"vvel") == 0) velocity.vvel = var;
    else if (strcmp(var_name,"wvel") == 0) velocity.wvel = var;
    else if (strcmp(var_name,"wgrd") == 0) velocity.wgrd = var;
    else if (strcmp(var_name,"btrc") == 0) velocity.btrc = var;
    else cerr << "unknown variable type = " << struct_name
              << "." << var_name << " undefined!" << endl;


  }

 

  else if (strcmp(struct_name,"temper") == 0) {
    if (strcmp(var_name,"temp") == 0) temper.temp = var;
    else if (strcmp(var_name,"bheatflx") == 0) temper.bheatflx = var;
    else if (strcmp(var_name,"bmlt") == 0) temper.bmlt = var;
    else cerr << "unknown variable type = " << struct_name
              << "." << var_name << " undefined!" << endl;
  }  
  else if (strcmp(struct_name,"numerics") == 0) { 
    if (strcmp(var_name,"tstart") == 0) numerics.tstart = var;    
    else if (strcmp(var_name,"tend") == 0) numerics.tend = var;    
    else if (strcmp(var_name,"time") == 0) numerics.time = var;    
    else cerr << "unknown variable type = " << struct_name
              << "." << var_name << " undefined!" << endl;
  }
  else if (strcmp(struct_name,"climate") == 0) {
    if (strcmp(var_name,"acab") == 0) climate.acab = var;
    else if (strcmp(var_name,"acab_tavg") == 0) climate.acab_tavg = var;
    else if (strcmp(var_name,"calving") == 0) climate.calving = var;
    else cerr << "unknown variable type = " << struct_name
              << "." << var_name << " undefined!" << endl;
  }
  else {
    cerr << "unknown variable type = " << struct_name
         << "." << var_name << " undefined!" << endl;
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
    else if (strcmp(var_name,"topg") == 0)
    {
      return(geometry.topg);
    }
    else if (strcmp(var_name,"usrf") == 0)
    {
      return(geometry.usrf);
    }
    else if (strcmp(var_name,"lsrf") == 0)
    {
      return(geometry.lsrf);
    }
    else if (strcmp(var_name,"floating_mask") == 0)
    {
      return(geometry.floating_mask);
    }
    else if (strcmp(var_name,"ice_mask") == 0)
    {
      return(geometry.ice_mask);
    }
    else if (strcmp(var_name,"lower_cell_loc") == 0)
    {
      return(geometry.lower_cell_loc);
    }
    else if (strcmp(var_name,"lower_cell_temp") == 0)
    {
      return(geometry.lower_cell_temp);
    }
    else  
      {
        cerr << "unknown variable type = " << struct_name
             << "." << var_name << " undefined!" << endl;
      }
  }
  else if (strcmp(struct_name,"numerics") == 0) {
    if (strcmp(var_name,"dew") == 0) return(numerics.dew);
    else if (strcmp(var_name,"dns") == 0) return(numerics.dns);
    else if (strcmp(var_name,"tstart") == 0) return(numerics.tstart);
    else if (strcmp(var_name,"tend") == 0) return(numerics.tend);
    else if (strcmp(var_name,"time") == 0) return(numerics.time);
    else {   
      cerr << "unknown variable type = " << struct_name
           << "." << var_name << " undefined!" << endl;
    }
  }
  else if (strcmp(struct_name,"constants") == 0) {
    if (strcmp(var_name,"seconds_per_year") == 0) return(&constants.seconds_per_year);
    else if (strcmp(var_name,"gravity") == 0) return(&constants.gravity);
    else if (strcmp(var_name,"rho_ice") == 0) return(&constants.rho_ice);
    else if (strcmp(var_name,"rho_seawater") == 0) return(&constants.rho_seawater);
    else if (strcmp(var_name,"therm_diffusivity_ice") == 0) return(&constants.therm_diffusivity_ice);
    else if (strcmp(var_name,"heat_capacity_ice") == 0) return(&constants.heat_capacity_ice);
    else {   
      cerr << "unknown variable type = " << struct_name
           << "." << var_name << " undefined!" << endl;
    }
  }
  else if (strcmp(struct_name,"velocity") == 0) {
    if (strcmp(var_name,"btrc") == 0) return (velocity.btrc);
    else if (strcmp(var_name,"uvel") == 0) return (velocity.uvel);
    else if (strcmp(var_name,"vvel") == 0) return (velocity.vvel);
    else if (strcmp(var_name,"wvel") == 0) return (velocity.wvel);
    else if (strcmp(var_name,"wgrd") == 0) return (velocity.wgrd);
    else {
      cerr << "unknown variable type = " << struct_name
           << "." << var_name << " undefined!" << endl;
    }
    //cout << "Set velocity var, " << var_name << endl;
  }

  else if (strcmp(struct_name,"temper") == 0) {
    if (strcmp(var_name,"temp") == 0) var = temper.temp;
    else if (strcmp(var_name,"bheatflx") == 0) var = temper.bheatflx;
    else if (strcmp(var_name,"bmlt") == 0) var = temper.bmlt;
    else {
      cerr << "unknown variable type = " << struct_name
           << "." << var_name << " undefined!" << endl;
    }
  }

  else if (strcmp(struct_name,"climate") == 0) {
    if (strcmp(var_name,"acab") == 0) var = climate.acab;
    else {
      cerr << "unknown variable type = " << struct_name
           << "." << var_name << " undefined!" << endl;
    }
  }
  else {
    cerr << "unknown variable type = " << struct_name
         << "." << var_name << " undefined!" << endl;
  }
  return(var);
}


int
DyCoreToGlimmer::setLongVar(long * var,  const char *var_name,  const char *struct_name)
{
  if (strcmp(struct_name,"geometry") == 0) {
    if (strcmp(var_name,"dimInfo") == 0) geometry.dimInfo = var;
  }
  else if (strcmp(struct_name,"velocity") == 0) {
    //cout << "Set velocity var, " << var_name << endl;
  }
  else {
    cerr << "unknown longVar type = " << struct_name
         << "." << var_name << endl;
  }
  return(0);
}



long * 
DyCoreToGlimmer::getLongVar( const char *var_name,  const char *struct_name)
{
  long * var;

  if (strcmp(struct_name,"geometry") == 0) {
    if (strcmp(var_name,"dimInfo") == 0) var = geometry.dimInfo;
    else if (strcmp(var_name,"ewlb") == 0) var = geometry.ewlb;
    else if (strcmp(var_name,"ewub") == 0) var = geometry.ewub;
    else if (strcmp(var_name,"nslb") == 0) var = geometry.nslb;
    else if (strcmp(var_name,"nsub") == 0) var = geometry.nsub;
    else if (strcmp(var_name,"nhalo") == 0) var = geometry.nhalo;
    else
      {
        cerr << "unknonwn variable " << var_name << " in "
             << struct_name << endl;
      }        
  }
  else if (strcmp(struct_name,"mpi_vars") == 0) {
    if (strcmp(var_name,"communicator") == 0) var = mpi_vars.communicator;
    else if (strcmp(var_name,"process_count") == 0) var = mpi_vars.process_count;
    else if (strcmp(var_name,"my_rank") == 0) var = mpi_vars.my_rank;
    else
      {
        cerr << "unknonwn variable " << var_name << " in "
             << struct_name << endl;
      }
  }

  else if (strcmp(struct_name,"velocity") == 0) {
    if (strcmp(var_name,"dimInfo") == 0) var = velocity.dimInfo;
    else
      {
        cerr << "unknonwn variable " << var_name << " in "
             << struct_name << endl;
      }
  }
  
  else if (strcmp(struct_name,"climate") == 0) {
    if (strcmp(var_name,"dimInfo") == 0) var = climate.dimInfo;
    else
      {
        cerr << "unknonwn variable " << var_name << " in "
             << struct_name << endl;
      }
  }
  
  return(var);
}


int
DyCoreToGlimmer::setInt4Var(int * var,  const char *var_name,  const char *struct_name)
{
  // cout << "struct_name::" << struct_name << "::" << endl;

  if (strcmp(struct_name,"felix_struct_name") == 0) {
    // if (strcmp(var_name,"dimInfo") == 0) geometry.dimInfo = var;
  }
  else if (strcmp(struct_name,"velocity") == 0) {
    //cout << "Set velocity var, " << var_name << endl;
  }
  else {
    cerr << "unknown int4Var type = " << struct_name
         << "." << var_name << endl;
  }
  return(0);
}


int * 
DyCoreToGlimmer::getInt4Var( const char *var_name,  const char *struct_name)
{
  int * var;

  if (strcmp(struct_name,"felix_struct_name") == 0) {
    if (strcmp(var_name,"dimInfo") == 0) var = 0;
    else
      {
        cerr << "unknonwn variable " << var_name << " in "
             << struct_name << endl;
      }        
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
    if (strcmp(var_name,"dew") == 0) 
      {
        numerics.dew = new double[elem_count];
        for (i=0;i<elem_count;i++) numerics.dew[i] = var[i];
        //cout << "Copy dew. dew, elem_count = "  << *numerics.dew << "  " << elem_count << endl;
      }  
    else if (strcmp(var_name,"dns") == 0) 
      {
        numerics.dns = new double[elem_count];
        for (i=0;i<elem_count;i++) numerics.dns[i] = var[i];
      //cout << "Copy dns. dns, elem_count = "  << *numerics.dns << "  " << elem_count << endl;
    }
    else 
      {
        cerr << "Unknown double variable name " << endl;
      }
  }
  
  // note that constants aren't pointers
  else if (strcmp(struct_name,"constants") == 0) {
    if (strcmp(var_name,"seconds_per_year") == 0) {
      constants.seconds_per_year = var[0];
    }
    else if (strcmp(var_name,"gravity") == 0) {
      constants.gravity = var[0];
    }
    else if (strcmp(var_name,"rho_ice") == 0) {
      constants.rho_ice = var[0];
    }
    else if (strcmp(var_name,"rho_seawater") == 0) {
      constants.rho_seawater = var[0];
    }
    else if (strcmp(var_name,"therm_diffusivity_ice") == 0) {
      constants.therm_diffusivity_ice = var[0];
    }
    else if (strcmp(var_name,"heat_capacity_ice") == 0) {
      constants.heat_capacity_ice = var[0];
    }    
    else
      {
        cerr << "Unknown double type: " << struct_name << "." 
             << var_name << endl;
      }
  }
  
  else if (strcmp(struct_name,"climate") == 0) {
    if (strcmp(var_name,"eus") == 0) {
      climate.eus = new double[elem_count];
      for (i=0;i<elem_count;i++) climate.eus[i] = var[i];
    }
    else
      {
        cerr << "Unknown climate double variable name: " << var_name << endl;
      }
  }
  else
    {
      cerr << "Unknown double structure name: " << struct_name << endl;
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
    else
      {
        cerr << "Unknown velocity integer var name: " << var_name << endl;
      }
  }
  else if (strcmp(struct_name,"geometry") == 0) {
    if (strcmp(var_name,"dimInfo") == 0) {
      geometry.dimInfo = new long[elem_count];
      for (i=0;i<elem_count;i++) geometry.dimInfo[i] = var[i];
    }  
    else if (strcmp(var_name,"ewlb") == 0) {
      geometry.ewlb = new long[elem_count];
      for (i=0;i<elem_count;i++) geometry.ewlb[i] = var[i];
    }  
    else if (strcmp(var_name,"ewub") == 0) {
      geometry.ewub = new long[elem_count];
      for (i=0;i<elem_count;i++) geometry.ewub[i] = var[i];
    }  
    else if (strcmp(var_name,"nslb") == 0) {
      geometry.nslb = new long[elem_count];
      for (i=0;i<elem_count;i++) geometry.nslb[i] = var[i];
    }  
    else if (strcmp(var_name,"nsub") == 0) {
      geometry.nsub = new long[elem_count];
      for (i=0;i<elem_count;i++) geometry.nsub[i] = var[i];
    }  
    else if (strcmp(var_name,"nhalo") == 0) {
      geometry.nhalo = new long[elem_count];
      for (i=0;i<elem_count;i++) geometry.nhalo[i] = var[i];
    }  
    else
      {
        cerr << "Unknown geometry integer variable name: " << var_name << endl;
      }
  }
  else if (strcmp(struct_name,"climate") == 0) {
    if (strcmp(var_name,"dimInfo") == 0) {
      climate.dimInfo = new long[elem_count];
      for (i=0;i<elem_count;i++) climate.dimInfo[i] = var[i];
    }  
    else
      {
        cerr << "Unknown climate integer variable name: " << var_name << endl;
      }
  }
  else if (strcmp(struct_name,"temper") == 0) {
    if (strcmp(var_name,"dimInfo") == 0) {
      temper.dimInfo = new long[elem_count];
      for (i=0;i<elem_count;i++) temper.dimInfo[i] = var[i];
    }  
    else
      {
        cerr << "Unknown temper integer variable name: " << var_name << endl;
      }
  }
  else if (strcmp(struct_name,"mpi_vars") == 0) {
    if (strcmp(var_name,"communicator") == 0) {
      mpi_vars.communicator = new long[elem_count];
      for (i=0;i<elem_count;i++) mpi_vars.communicator[i] = var[i];
    }
    else if (strcmp(var_name,"process_count") == 0) {
      mpi_vars.process_count = new long[elem_count];
      for (i=0;i<elem_count;i++) mpi_vars.process_count[i] = var[i];
    }
    else if (strcmp(var_name,"my_rank") == 0) {
      mpi_vars.my_rank = new long[elem_count];
      for (i=0;i<elem_count;i++) mpi_vars.my_rank[i] = var[i];
    }
    else
      {
        cerr << "Unknown mpi_vars integer variable name: " << var_name << endl;
      }
  }
  else 
    {
      cerr << "Unknown integer struc_name: " << struct_name << endl;
    }
  
  return(0);
}


int 
DyCoreToGlimmer::initDyCore(const char * dycore_fname)
{
  // cout << "In DycoreToGlimmer::initDyCore" << endl;
  return(0);
}

int
DyCoreToGlimmer::runDyCore(double& cur_time_yr, const double time_inc_yr)
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
