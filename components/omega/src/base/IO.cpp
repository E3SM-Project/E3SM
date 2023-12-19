//===-- base/IO.cpp - basic IO utilities implementation ---------*- C++ -*-===//
//
// The routines here provide the interfaces to the parallel IO environment,
// currently the SCORPIO parallel IO library.
// These functions should generally be accessed through the IOStreams
// interfaces and not used directly.
//
//===----------------------------------------------------------------------===//

#include "IO.h"
#include "DataTypes.h"
#include "Logging.h"
#include "mpi.h"
#include "pio.h"

#include <map>
#include <string>

namespace OMEGA {

// Define global variables
//------------------------------------------------------------------------------
int IOSysID                = 0;
IOFileFmt IODefaultFileFmt = IOFmtDefault;

// Utilities
//------------------------------------------------------------------------------
// Converts string choice for PIO rearranger to an enum
IORearranger
IORearrFromString(const std::string &Rearr // [in] choice of IO rearranger
) {
   // Set default return value
   IORearranger ReturnRearr = IORearrUnknown;

   // Convert input string to lowercase for easier comparison
   std::string RearrComp = Rearr;
   std::transform(RearrComp.begin(), RearrComp.end(), RearrComp.begin(),
                  [](unsigned char c) { return std::tolower(c); });

   // Determine the appropriate enum to use based on the name
   // Check most likely options first

   // Box rearranger
   if (RearrComp == "box") {
      ReturnRearr = IORearrBox;

      // Subset rearranger
   } else if (RearrComp == "subset") {
      ReturnRearr = IORearrSubset;

      // Default
   } else if (RearrComp == "default") {
      ReturnRearr = IORearrDefault;

   } else {
      ReturnRearr = IORearrUnknown;
   }

   return ReturnRearr;

} // End IORearrFromString

//------------------------------------------------------------------------------
// Converts string choice for File Format to an enum
IOFileFmt
IOFileFmtFromString(const std::string &Format // [in] choice of IO file format
) {
   // Set default return value
   IOFileFmt ReturnFileFmt = IOFmtUnknown;

   // Convert input string to lowercase for easier comparison
   std::string FmtCompare = Format;
   std::transform(FmtCompare.begin(), FmtCompare.end(), FmtCompare.begin(),
                  [](unsigned char c) { return std::tolower(c); });

   // Determine the appropriate enum to use based on the input string
   // Check most likely options first

   // NetCDF4 variants should be default so check first
   // we assume netcdf4 refers to netcdf4p
   if (FmtCompare == "netcdf4") {
      ReturnFileFmt = IOFmtNetCDF4;

      // Check for ADIOS
   } else if (FmtCompare == "adios") {
      ReturnFileFmt = IOFmtADIOS;

      // Check specific netcdf4 variants - compressed
   } else if (FmtCompare == "netcdf4c") {
      ReturnFileFmt = IOFmtNetCDF4c;

      // Check specific netcdf4 variants - parallel
   } else if (FmtCompare == "netcdf4p") {
      ReturnFileFmt = IOFmtNetCDF4p;

      // Check for older netcdf
   } else if (FmtCompare == "netcdf3") {
      ReturnFileFmt = IOFmtNetCDF3;

      // Check for older pnetcdf
   } else if (FmtCompare == "pnetcdf") {
      ReturnFileFmt = IOFmtPnetCDF;

      // Check for native HDF5
   } else if (FmtCompare == "hdf5") {
      ReturnFileFmt = IOFmtHDF5;

   } else if (FmtCompare == "default") {
      ReturnFileFmt = IOFmtDefault;

   } else {
      ReturnFileFmt = IOFmtUnknown;
   }

   return ReturnFileFmt;

} // end of IOFileFmtFromString

//------------------------------------------------------------------------------
// Converts string choice for IO read/write mode to an enum
IOMode
IOModeFromString(const std::string &Mode // [in] choice of IO mode (read/write)
) {
   // Set default return value
   IOMode ReturnMode = IOModeUnknown;

   // Convert input string to lowercase for easier comparison
   std::string ModeComp = Mode;
   std::transform(ModeComp.begin(), ModeComp.end(), ModeComp.begin(),
                  [](unsigned char c) { return std::tolower(c); });

   // Determine the appropriate enum to use based on the input string

   // Read only
   if (ModeComp == "read") {
      ReturnMode = IOModeRead;

      // Write only
   } else if (ModeComp == "write") {
      ReturnMode = IOModeWrite;

   } else {
      ReturnMode = IOModeUnknown;
   }

   return ReturnMode;

} // End IOModeFromString

//------------------------------------------------------------------------------
// Converts string choice for existence behavior to an enum
IOIfExists IOIfExistsFromString(
    const std::string &IfExists // [in] choice of behavior on file existence
) {
   // Set default return value
   IOIfExists ReturnIfExists = IOIfExists::Fail;

   // Convert input string to lowercase for easier comparison
   std::string IECompare = IfExists;
   std::transform(IECompare.begin(), IECompare.end(), IECompare.begin(),
                  [](unsigned char c) { return std::tolower(c); });

   // Determine the appropriate enum to use based on the input string

   // Fail with error
   if (IECompare == "fail") {
      ReturnIfExists = IOIfExists::Fail;

      // Replace existing file
   } else if (IECompare == "replace") {
      ReturnIfExists = IOIfExists::Replace;

      // Append or insert into existing file
   } else if (IECompare == "append") {
      ReturnIfExists = IOIfExists::Append;

   } else {
      ReturnIfExists = IOIfExists::Fail;
   }

   return ReturnIfExists;

} // End IOIfExistsFromString

//------------------------------------------------------------------------------
// Converts string choice for floating point precision to an enum
IOPrecision IOPrecisionFromString(
    const std::string &Precision // [in] choice of floating point precision
) {
   // Set default return value
   IOPrecision ReturnPrec = IOPrecision::Double;

   // Convert input string to lowercase for easier comparison
   std::string PrecComp = Precision;
   std::transform(PrecComp.begin(), PrecComp.end(), PrecComp.begin(),
                  [](unsigned char c) { return std::tolower(c); });

   // Determine the appropriate enum to use based on the string name

   // Maintain full double precision
   if (PrecComp == "double") {
      ReturnPrec = IOPrecision::Double;

      // Write reduced (single) precision
   } else if (PrecComp == "single") {
      ReturnPrec = IOPrecision::Single;

   } else {
      ReturnPrec = IOPrecision::Double;
   }

   return ReturnPrec;

} // End IOPrecisionFromString

// Methods
//------------------------------------------------------------------------------
// Initializes the IO system based on configuration inputs and
// default MPI communicator
int IOInit(const MPI_Comm &InComm // [in] MPI communicator to use
) {

   int Err = 0; // success error code
   // Retrieve parallel IO parameters from the Omega configuration
   // TODO: replace with actual config, for now hardwired
   // extern int IOSysID;

   IOFileFmt IODefaultFileFmt = IOFileFmtFromString("netcdf4c");
   int NumIOTasks             = 1;
   int IOStride               = 1;
   int IOBaseTask             = 0;
   IORearranger Rearrange     = IORearrDefault;

   // Call PIO routine to initialize
   Err = PIOc_Init_Intracomm(InComm, NumIOTasks, IOStride, IOBaseTask,
                             Rearrange, &IOSysID);
   if (Err != 0)
      LOG_ERROR("IOInit: Error initializing SCORPIO");

   return Err;

} // end IOInit

//------------------------------------------------------------------------------
// This routine opens a file for reading or writing, depending on the
// Mode argument. The filename with full path must be supplied and
// a FileID is returned to be used by other IO functions.
// The format of the file is assumed to be the default defined on init
// but can be optionally changed through this open function.
// For files to be written, optional arguments govern the behavior to be
// used if the file already exists, and the precision of any floating point
// variables. Returns an error code.
int IOFileOpen(
    int &FileID,                 // [out] returned fileID for this file
    const std::string &Filename, // [in] name (incl path) of file to open
    IOMode Mode,                 // [in] mode (read or write)
    IOFileFmt InFormat,          // [in] (optional) file format
    IOIfExists IfExists,         // [in] (for writes) behavior if file exists
    IOPrecision Precision        // [in] (for writes) precision of floats
) {

   int Err    = 0;        // default success return code
   int Format = InFormat; // coerce to integer for PIO calls

   switch (Mode) {

   // If reading, open the file for read-only
   case IOModeRead:
      Err = PIOc_openfile(IOSysID, &FileID, &Format, Filename.c_str(), Mode);
      if (Err != PIO_NOERR)
         LOG_ERROR("IOFileOpen: PIO error opening file {} for read", Filename);
      break;

   // If writing, the open will depend on what to do if the file
   // exists.
   case IOModeWrite:

      switch (IfExists) {
      // If the write should be a new file and fail if the
      // file exists, we use create and fail with an error
      case IOIfExists::Fail:
         Err = PIOc_createfile(IOSysID, &FileID, &Format, Filename.c_str(),
                               NC_NOCLOBBER | Mode);
         if (Err != PIO_NOERR)
            LOG_ERROR("IOFileOpen: PIO error opening file {} for writing",
                      Filename);
         break;

      // If the write should replace any existing file
      // we use create with the CLOBBER option
      case IOIfExists::Replace:
         Err = PIOc_createfile(IOSysID, &FileID, &Format, Filename.c_str(),
                               NC_CLOBBER | Mode);
         if (Err != PIO_NOERR)
            LOG_ERROR("IOFileOpen: PIO error opening file {} for writing",
                      Filename);
         break;

      // If the write should append or add to an existing file
      // we open the file for writing
      case IOIfExists::Append:
         Err = PIOc_openfile(IOSysID, &FileID, &Format, Filename.c_str(), Mode);
         if (Err != PIO_NOERR)
            LOG_ERROR("IOFileOpen: PIO error opening file {} for writing",
                      Filename);
         break;

      default:
         LOG_ERROR("IOFileOpen: unknown IfExists option for writing");

      } // end switch IfExists
      break;

   // Unknown mode
   default:
      LOG_ERROR("IOFileOpen: Unknown Mode for file");

   } // End switch on Mode

   return Err;

} // End IOFileOpen
//------------------------------------------------------------------------------
// Closes an open file using the fileID, returns an error code
int IOFileClose(int &FileID /// [in] ID of the file to be closed
) {
   // Just calls the PIO close routine
   int Err = PIOc_closefile(FileID);
   return Err;

} // End IOFileClose

//------------------------------------------------------------------------------
// Retrieves a dimension length from an input file, given the name
// of the dimension. Returns the length if exists, but returns a negative
// value if a dimension of that length is not found in the file.
int IOGetDimLength(int FileID, // [in] ID of the file containing dim
                   const std::string &DimName // [in] existing parent MachEnv
) {

   int Err = 0; // Local error code

   // First get the dimension ID
   int DimID = 0;
   Err       = PIOc_inq_dimid(FileID, DimName.c_str(), &DimID);
   if (Err != PIO_NOERR) {
      LOG_ERROR("IOGetDimLength: PIO error while retrieving dimensionID");
      return -1;
   }

   // Now retrieve the length and return
   PIO_Offset Length = -1;
   Err               = PIOc_inq_dimlen(FileID, DimID, &Length);
   if (Err != PIO_NOERR) {
      LOG_ERROR("IOGetDimLength: PIO error while retrieving dimension length");
      return -1;
   }
   return Length;

} // End IOGetDimLength

//------------------------------------------------------------------------------
// Creates a PIO decomposition description to describe the layout of
// a distributed array of given type.
int IOCreateDecomp(
    int &DecompID,      // [out] ID assigned to the new decomposition
    IODataType VarType, // [in] data type of array
    int NDims,          // [in] number of array dimensions
    const std::vector<int> &DimLengths, // [in] global dimension lengths
    int Size,                           // [in] local size of array
    const std::vector<int> &GlobalIndx, // [in] global indx for each local indx
    IORearranger Rearr                  // [in] rearranger method to use
) {

   int Err = 0; // default return code

   // Convert global index array into an offset array expected by PIO
   std::vector<PIO_Offset> CompMap;
   CompMap.resize(Size);
   for (int i = 0; i < Size; ++i) {
      CompMap[i] = GlobalIndx[i];
   } // end global index loop

   // Call the PIO routine to define the decomposition
   // int TmpRearr = Rearr; // needed for type compliance across interface
   // Err = PIOc_InitDecomp(IOSysID, VarType, NDims, DimLengths, Size, CompMap,
   //                      &DecompID, &TmpRearr, nullptr, nullptr);
   Err = PIOc_init_decomp(IOSysID, VarType, NDims, &DimLengths[0], Size,
                          &CompMap[0], &DecompID, Rearr, nullptr, nullptr);
   if (Err != PIO_NOERR)
      LOG_ERROR("IOCreateDecomp: PIO error defining decomposition");

   return Err;

} // End IOCreateDecomp

//------------------------------------------------------------------------------
// Removes a defined PIO decomposition description to free memory
int IODestroyDecomp(int &DecompID // [inout] ID for decomposition to be removed
) {

   int Err = PIOc_freedecomp(IOSysID, DecompID);
   if (Err != PIO_NOERR)
      LOG_ERROR("IODestroyDecomp: PIO error freeing decomposition");

   return Err;

} // End IODestroyDecomp

//------------------------------------------------------------------------------
// Reads a distributed array. This overloaded interface is for integer
// arrays. Also, all arrays are assumed to be in contiguous storage so
// the arrays of any dimension are treated as a 1-d array with the full
// local size.
int IOReadArray(int *Array,                 ///< [out] array to be read
                int Size,                   ///< [in] local size of array
                const std::string &VarName, ///< [in] name of variable to read
                int FileID,  ///< [in] ID of open file to read from
                int DecompID ///< [in] decomposition ID for this var
) {

   int Err   = 0; // default return code
   int VarID = 0;

   // Find variable ID from file
   Err = PIOc_inq_varid(FileID, VarName.c_str(), &VarID);
   if (Err != PIO_NOERR)
      LOG_ERROR("IOReadArray(int): Error finding varid");

   // PIO Read array call to read the distributed array
   void *ArrayPtr   = Array;
   PIO_Offset ASize = Size;
   Err = PIOc_read_darray(FileID, VarID, DecompID, ASize, ArrayPtr);
   if (Err != PIO_NOERR)
      LOG_ERROR("IOReadArray(int): Error in SCORPIO read array");

   return Err;

} // End IOReadArray (int)

//------------------------------------------------------------------------------
// Reads a distributed array. This overloaded interface is for 32-bit real
// arrays. Also, all arrays are assumed to be in contiguous storage so
// the arrays of any dimension are treated as a 1-d array with the full
// local size.
int IOReadArray(float *Array,               ///< [out] array to be read
                int Size,                   ///< [in] local size of array
                const std::string &VarName, ///< [in] name of variable to read
                int FileID,  ///< [in] ID of open file to read from
                int DecompID ///< [in] decomposition ID for this var
) {

   int Err   = 0; // default return code
   int VarID = 0;

   // Find variable ID from file
   Err = PIOc_inq_varid(FileID, VarName.c_str(), &VarID);
   if (Err != PIO_NOERR)
      LOG_ERROR("IOReadArray(real): Error finding varid");

   // PIO Read array call to read the distributed array
   Err = PIOc_read_darray(FileID, VarID, DecompID, Size, Array);
   if (Err != PIO_NOERR)
      LOG_ERROR("IOReadArray(real): Error in SCORPIO read array");

   return Err;

} // End IOReadArray (real)

//------------------------------------------------------------------------------
// Reads a distributed array. This overloaded interface is for 64-bit real
// arrays. Also, all arrays are assumed to be in contiguous storage so
// the arrays of any dimension are treated as a 1-d array with the full
// local size.
int IOReadArray(double *Array,              ///< [out] array to be read
                int Size,                   ///< [in] local size of array
                const std::string &VarName, ///< [in] name of variable to read
                int FileID,  ///< [in] ID of open file to read from
                int DecompID ///< [in] decomposition ID for this var
) {

   int Err   = 0; // default return code
   int VarID = 0;

   // Find variable ID from file
   Err = PIOc_inq_varid(FileID, VarName.c_str(), &VarID);
   if (Err != PIO_NOERR)
      LOG_ERROR("IOReadArray(double): Error finding varid");

   // PIO Read array call to read the distributed array
   Err = PIOc_read_darray(FileID, VarID, DecompID, Size, Array);
   if (Err != PIO_NOERR)
      LOG_ERROR("IOReadArray(double): Error in SCORPIO read array");

   return Err;

} // End IOReadArray (double)

//------------------------------------------------------------------------------

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
