// The DyCoreToGlimmer class provides methods to move Glimmer Fortran data to C++ structures
// for access by C++ based dynamical cores.  The structure names and structure member
// names mostly correspond to derived types defined in Glimmer.  In general, pointers to
// the Fortran data arrays are used, rather than copies of these arrays.  This saves space
// and reduces the steps needed to update the Glimmer data between calls to the core
// ice sheet modeling program.  Methods are provided to set these array pointers, and 
// copy array dimension information.  Objects of this class are accessed by extern C
// routines in dycore_to_glimmer_extern.cpp, and by the dynamical core front end.  DMR--5/24/10

//#pragma once
#ifndef DYCORETOGLIMMER
#define DYCORETOGLIMMER

#include <stdio.h>
#include <iostream>


class DyCoreToGlimmer
{

  private:

  // Keep track of dycore_type and dycore_index.  Dycore_index is used to index the
  // the external dycore object storage array.    
  struct {
    int dycore_type;
    int dycore_index;
  } dycore_info;

  // The following structures are based on the derived types in glide_types.F90 

  struct {
    double * thck;
    double * usrf;
    double * lsrf;
    double * topg;
    double * floating_mask;
    double * ice_mask;
    double * lower_cell_loc; // z-location of lowest cell-center
    double * lower_cell_temp; // temperature in lowest cell 
    long * dimInfo;
    long * ewlb;
    long * ewub;
    long * nslb;
    long * nsub;
    long * nhalo;

    //double * thkmask;
    //double * marine_bc_normal;
 
  } geometry;

  struct {
    double * uvel; //output
    double * vvel; //output
    double * wvel;
    double * wgrd;
    double * btrc; // basal traction coefficient
    long * dimInfo;
  } velocity;

  struct {
    double * temp; // Three-dimensional temperature field.
    double * bheatflx; // basal heat flux (2D)
    double * bmlt;  // Basal melt-rate
    long * dimInfo;
  } temper;

  struct {

  } lithot_type;

  struct {
    double * tstart;
    double * tend;
    double * time; 
    //double * tinc;  

    double * dew; // ew cell size
    double * dns; // ns cell size
    
  } numerics;

  struct {
    double * acab; // Annual mass balance.
    double * acab_tavg; // Annual mass balance (time average)
    double * calving; // Calving flux (scaled as mass balance, thickness, etc)
    long * dimInfo;
    double * eus; // eustatic sea level
  } climate;

  struct { 
    double * beta;   // basal shear coefficient
    double * btraction; // -dir (1,:,:) and y-dir (2,:,:) "consistent" basal 
                        // traction fields (calculated from matrix coeffs)
    long dimInfo;
  } velocity_hom;

  struct {
    double seconds_per_year;
    double gravity;
    double rho_ice; 
    double rho_seawater;  
    double therm_diffusivity_ice;
    double heat_capacity_ice;
  } constants;


  struct {
    long * communicator;
    long * process_count;
    long * my_rank;
  } mpi_vars;

  public: 

  DyCoreToGlimmer();
  virtual ~DyCoreToGlimmer();

  
  int setDoubleVar( double *var, const char *var_name,  const char *struct_name);
  double * getDoubleVar( const char *var_name,  const char *struct_name);

  int setLongVar( long * var,  const char * var_name,  const char *struct_name);
  long * getLongVar( const char * var_name,  const char *struct_name);

  int setInt4Var( int * var,  const char * var_name,  const char *struct_name);
  int * getInt4Var( const char * var_name,  const char *struct_name);

  int copyInDoubleVar( const double *var,  const char *var_name, 
		        const char *struct_name,  const long *var_dim_info);
  int copyInLongVar( const long *var,  const char *var_name, 
		      const char *struct_name,  const long *var_dim_info);
    
  virtual int initDyCore(const char *input_fname);    // = 0;
  virtual int runDyCore(double& cur_time_yr, const double time_inc_yr);     // = 0;
  virtual int deleteDyCore(); // = 0; 

  int setDyCoreType(const int dycore_type);
  int getDyCoreType();

  int setDyCoreIndex(const int dycore_index);
  int getDyCoreIndex();

};

#endif
