// $Id: compare.cpp 5181 2009-07-14 11:09:04Z martin-johnson $

// Program to compare any variables common to two netcdf files.
// Optional tolerances:-
// Absolute:
// if the difference between values A and B is less than the
// absolute tolerance, then it is small enough to ignore.
// Relative:
// If the difference between values A and B is less than or
// equal to the given number of Ulps, then they are close
// enough.

#include <stdlib.h>
#include <iostream>
#include <cstdio>
#include <cstring>
#include <cfloat>
#include <cstring>
#include <cmath>
#include <cstring>
#include <cassert>
#include <netcdfcpp.h>

#if !defined(WIN32)
#include <unistd.h>
#else
// unistd is a UNIX specific library so we must provide an alternate getopt() function
// Prototype the getopt function provided in the NetCDF distribution (src/win32/NET)
extern int optind;
extern int getopt(int argc,  char *const argv[], const char *opts);
extern char *optarg;
#endif

using namespace std;

// option parsing
// -v for verbose
// -s for silent
// -z to fail if either file contains no variables
#define OPTIONS "vszr:a:"
bool verbose = false;
bool silent = false;
bool zerocheck = false;

// Define some error types
#define ERRDIMNUM 1
#define ERRTYPE   2
#define ERRSIZE   3
#define ERRVALS   4

// Prototype comparison functions
int varComp(NcVar *A, NcVar *B, double dMinAbs, int iMaxUlps);
bool AlmostEqual2sComplement(float A, float B, int maxUlps);
int twosComplementDiff(float A, float B);

int main(int argc, char **argv)
{

  int iMaxUlps = 0;          // max num of Ulps allowed between floats A and B
  double dMinAbs = DBL_MIN;  // min absolute difference between floats allowed

  NcToken* namesAP;
  NcToken* namesBP;
  int namesCommon;

  bool bErrFound = false;  // flag if at least one comparison error has been found

  // Parse options
  char c;
  while ((c = getopt(argc, argv, OPTIONS)) != -1)
    switch (c) {
	  case 'v':
	      verbose=true;
	      break;
	  case 's':
	      silent=true;
	      break;
	  case 'z':
	      zerocheck=true;
	      break;
          case 'r':
              iMaxUlps = atoi(optarg);
              break;
          case 'a':
	      dMinAbs = atof(optarg);
              break;
      }

  if (silent) verbose=false;

  // check we have some options left for filenames and tolerance
  if (argc-optind != 2) {
    cout << "Usage: compare [options] <NCfileA> <NCfileB>" << endl;
    cout << "Options:" << endl;
    cout << "  -v\t\t\tVerbose mode" << endl;
    cout << "  -s\t\t\tSilent mode" << endl;
    cout << "  -z\t\t\tFail for empty files" << endl;
    cout << "  -r <val>\t\tRelative tolerance threshold in ulps" << endl;
    cout << "  -a <val>\t\tAbsolute tolerance threshold--ignore smaller numbers" << endl;
    exit(1);
  }

  const char *fileA = argv[optind];
  const char *fileB = argv[1+optind];

  // ---------------------------------------------------------------
  // Verify the threshold options, if supplied 
  // ---------------------------------------------------------------

  if(iMaxUlps<0) {
    fprintf(stderr,"Error: Cannot specify a -ve relative tolerance\n");
    exit(EXIT_FAILURE);
  }
  if(dMinAbs<0.0) {
    fprintf(stderr,"Error: Cannot specify a -ve absolute tolerance\n");
    exit(EXIT_FAILURE);
  }

  // error handling class
  NcError err_handler;

  // ---------------------------------------------------------------
  // Load the two files and their variable names 
  // ---------------------------------------------------------------

  // == NC file A ==
  // Check existance
  NcFile fpA(fileA);
  if(!fpA.is_valid()) {
    fprintf(stderr,"Error: Could not open file: %s, or it is not a valid NetCDF file\n", fileA);
    exit(EXIT_FAILURE);
  }
  // get number of variables
  const int numVarsA = fpA.num_vars();

  // check there are variables
  if (zerocheck && (numVarsA==0)) {
      fprintf(stderr,"File %s contains no variables\n",fileA);
    exit(EXIT_FAILURE);
  }

  // allocate space to store variable names
  namesAP = (NcToken*)malloc(sizeof(NcToken)*numVarsA);

  // get variable names
  for(int ii=0; ii<numVarsA; ++ii) {
    namesAP[ii] = fpA.get_var(ii)->name();
  }

  // == NC file B ==
  // Check existance
  NcFile fpB(fileB);
  if(!fpB.is_valid()) {
    fprintf(stderr,"Error: Could not open file: %s, or it is not a valid NetCDF file\n", fileB);
    exit(EXIT_FAILURE);
  }
  // get number of variables
  const int numVarsB = fpB.num_vars();

  // check there are variables
  if (zerocheck && (numVarsB==0)) {
      fprintf(stderr,"File %s contains no variables\n",fileB);
    exit(EXIT_FAILURE);
  }

  // allocate space to store variable names
  namesBP = (NcToken*)malloc(sizeof(NcToken)*numVarsB);

  // get variable names
  for(int ii=0; ii<numVarsB; ++ii) {
    namesBP[ii] = fpB.get_var(ii)->name();
  }

  // Change the error behavior of the netCDF C++ API by creating an
  // NcError object. Until it is destroyed, this NcError object will
  // ensure that the netCDF C++ API silently returns error codes on
  // any failure, and leaves any other error handling to the calling
  // program. In the case of this example, we just exit with an
  // NC_ERR error code.
  NcError err(NcError::silent_nonfatal);

  // -------------------------------------------------------------
  // Determine which variables in A are also
  // present in B, then for those that are, compare them
  // -------------------------------------------------------------

  for(int ii=0; ii<numVarsA; ++ii) {
      namesCommon=0;
      // Identify matching names in file B
      for (int jj=0; jj<numVarsB; ++jj) {
	  if (!strcmp(namesAP[ii],namesBP[jj])) {
	      namesCommon=jj;
	      break;
	  }
      }

      // Compare matching variables
      if (namesCommon) {
	  if (verbose)
	      cout << "Comparing variable " << namesAP[ii] << endl;
	  switch (varComp(fpA.get_var(ii),fpB.get_var(namesCommon),dMinAbs,iMaxUlps)) { 
	      case ERRDIMNUM:
		bErrFound = true;
		  if (verbose)
		      fprintf(stderr,"**ERROR: Differing numbers of dimensions in variable %s\n",namesAP[ii]);
		  break;
	      case ERRTYPE:
		bErrFound = true;
		  if (verbose)
		      fprintf(stderr,"**ERROR: Differing types in variable %s\n",namesAP[ii]);
		  break;
	      case ERRSIZE:
		bErrFound = true;
		  if (verbose)
		      fprintf(stderr,"**ERROR: Differing sizes of variable %s\n",namesAP[ii]);
		  break;
	      case ERRVALS:
		bErrFound = true;
		  if (verbose)
		      fprintf(stderr,"**ERRROR: Differing values in variable %s\n",namesAP[ii]);
		  break;
	      default:
		  continue;
	  }
      }
  }
  if (bErrFound) {
    if (!silent)
      fprintf(stderr,"Files %s and %s differ\n", fileA, fileB);
    exit(EXIT_FAILURE);
  }

  // free up allocated space
  free(namesAP);
  free(namesBP);

  return EXIT_SUCCESS;
}

// The comparison function -- returns 0 if identical, an error code if they differ
int varComp(NcVar *A, NcVar *B, double dMinAbs, int iMaxUlps)
{
    // First check types and dimensions
    if (A->num_dims()!=B->num_dims()) 
	return ERRDIMNUM;
    if (A->type()!=B->type()) 	
	return ERRTYPE;

    // Get all the values
    NcValues* AvalsP=A->values();
    NcValues* BvalsP=B->values();

    // Check sizes of variables
    long numValsA=AvalsP->num();
    long numValsB=BvalsP->num();

    if (numValsA!=numValsB) 
	return ERRSIZE;

    // Check if variables are scaled
    double ScaleA = 1.0;
    NcAtt* ScaleFactorA;
    if ((ScaleFactorA=A->get_att("scale_factor")))
      ScaleA = ScaleFactorA->as_double(0);
    double ScaleB = 1.0;
    NcAtt* ScaleFactorB;
    if ((ScaleFactorB=B->get_att("scale_factor")))
      ScaleB = ScaleFactorB->as_double(0);

    // Compare individual values
    int errcatch=0;
    for (int ii=0; ii<numValsA; ++ii) {
      // local copy of the values
      double Aval = ScaleA*AvalsP->as_double(ii);
      double Bval = ScaleB*BvalsP->as_double(ii);
      // if equal, then happy days
      if (Aval == Bval) {
	continue;
      }
      else {
	// if an absolute threshold is specified for ignoring
	if (dMinAbs > DBL_MIN) {
	  if (fabs(Aval - Bval) < dMinAbs) {
	    if (verbose) {
	      cout << Aval << " - " << Bval << " < " << dMinAbs << " (absolute threshold) and so skipping" << endl;
	    }
	    continue;
	  }
	}
	// if a relative threhold is specified for comparing
	if (iMaxUlps > 0) {
	  if (AlmostEqual2sComplement(Aval,Bval,iMaxUlps)) {
	    if (verbose) {
	      cout << Aval << " - " << Bval << " has an AlmostEqual2sComplement <= " << iMaxUlps 
		   << " (relative threshold) and so skipping" << endl;
	    }
	    continue;
	  }
	}
	// catch all is that Aval != Bval and so we have a problem
	errcatch=(errcatch || 1);
	if (verbose)
	  cout << Aval << " - " << Bval << " = " << Aval-Bval << "\t(" 
	       << twosComplementDiff(Aval,Bval) << " Ulps difference)" << endl;
	break;
      }
    }
    
    free(AvalsP);
    free(BvalsP);
    free(ScaleFactorA);
    free(ScaleFactorB);

    return errcatch*ERRVALS;
}

// Usable AlmostEqual function
bool AlmostEqual2sComplement(float A, float B, int maxUlps)
{
  // Make sure maxUlps is non-negative and small enough that the
  // default NAN won't compare as equal to anything.
  assert(maxUlps > 0 && maxUlps < 4 * 1024 * 1024);
  int aInt = *(int*)&A;

  // Make aInt lexicographically ordered as a twos-complement int
  if (aInt < 0)
    aInt = 0x80000000 - aInt;
  
  // Make bInt lexicographically ordered as a twos-complement int
  int bInt = *(int*)&B;
  if (bInt < 0)
    bInt = 0x80000000 - bInt;
  int intDiff = abs(aInt - bInt);
  if (intDiff <= maxUlps)
    return true;
  return false;
}

// Just the raw 2sComplement difference
int twosComplementDiff(float A, float B)
{
  int aInt = *(int*)&A;

  // Make aInt lexicographically ordered as a twos-complement int
  if (aInt < 0)
    aInt = 0x80000000 - aInt;
  
  // Make bInt lexicographically ordered as a twos-complement int
  int bInt = *(int*)&B;
  if (bInt < 0)
    bInt = 0x80000000 - bInt;
  int intDiff = abs(aInt - bInt);

  return intDiff;
}
