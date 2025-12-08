//===-- base/IO.cpp - basic IO utilities implementation ---------*- C++ -*-===//
//
// The routines here provide the interfaces to the parallel IO environment,
// currently the SCORPIO parallel IO library.
// These functions should generally be accessed through the IOStreams
// interfaces and not used directly.
//
//===----------------------------------------------------------------------===//

#include "IO.h"
#include "Config.h"
#include "DataTypes.h"
#include "Error.h"
#include "Logging.h"
#include "mpi.h"
#include "pio.h"

#include <filesystem>
#include <map>
#include <string>

namespace OMEGA {
namespace IO {

// Define global variables
//------------------------------------------------------------------------------
int SysID               = 0;
FileFmt DefaultFileFmt  = FmtDefault;
Rearranger DefaultRearr = RearrDefault;

// Utilities
//------------------------------------------------------------------------------
// Converts string choice for PIO rearranger to an enum
Rearranger
RearrFromString(const std::string &Rearr // [in] choice of IO rearranger
) {
   // Set default return value
   Rearranger ReturnRearr = RearrUnknown;

   // Convert input string to lowercase for easier comparison
   std::string RearrComp = Rearr;
   std::transform(RearrComp.begin(), RearrComp.end(), RearrComp.begin(),
                  [](unsigned char c) { return std::tolower(c); });

   // Determine the appropriate enum to use based on the name
   // Check most likely options first

   if (RearrComp == "box") { // Box rearranger
      ReturnRearr = RearrBox;

   } else if (RearrComp == "subset") { // Subset rearranger
      ReturnRearr = RearrSubset;

   } else if (RearrComp == "default") { // Default
      ReturnRearr = RearrDefault;

   } else {
      ABORT_ERROR("IO: Unknown data rearranger {}", RearrComp);
   }

   return ReturnRearr;

} // End RearrFromString

//------------------------------------------------------------------------------
// Converts string choice for File Format to an enum
FileFmt
FileFmtFromString(const std::string &Format // [in] choice of IO file format
) {
   // Set default return value
   FileFmt ReturnFileFmt = FmtUnknown;

   // Convert input string to lowercase for easier comparison
   std::string FmtCompare = Format;
   std::transform(FmtCompare.begin(), FmtCompare.end(), FmtCompare.begin(),
                  [](unsigned char c) { return std::tolower(c); });

   // Determine the appropriate enum to use based on the input string
   // Check most likely options first

   if (FmtCompare == "netcdf4") { // NetCDF4 variants
      ReturnFileFmt = FmtNetCDF4;

   } else if (FmtCompare == "adios") { // ADIOS
      ReturnFileFmt = FmtADIOS;

   } else if (FmtCompare == "netcdf4c") { // netcdf4 variant - compressed
      ReturnFileFmt = FmtNetCDF4c;

   } else if (FmtCompare == "netcdf4p") { // netcdf4 variant - parallel
      ReturnFileFmt = FmtNetCDF4p;

   } else if (FmtCompare == "netcdf3") { // Older netcdf
      ReturnFileFmt = FmtNetCDF3;

   } else if (FmtCompare == "pnetcdf") { // pnetcdf
      ReturnFileFmt = FmtPnetCDF;

   } else if (FmtCompare == "hdf5") { // native HDF5
      ReturnFileFmt = FmtHDF5;

   } else if (FmtCompare == "default") {
      ReturnFileFmt = FmtDefault;

   } else {
      ABORT_ERROR("IO: Unknown default file format {}", FmtCompare);
   }

   return ReturnFileFmt;

} // end of FileFmtFromString

//------------------------------------------------------------------------------
// Converts string choice for IO read/write mode to an enum
Mode ModeFromString(
    const std::string &ModeChoice // [in] IO mode choice (read/write)
) {
   // Set default return value
   Mode ReturnMode = ModeUnknown;

   // Convert input string to lowercase for easier comparison
   std::string ModeComp = ModeChoice;
   std::transform(ModeComp.begin(), ModeComp.end(), ModeComp.begin(),
                  [](unsigned char c) { return std::tolower(c); });

   // Determine the appropriate enum to use based on the input string

   if (ModeComp == "read") { // Read only
      ReturnMode = ModeRead;

   } else if (ModeComp == "write") { // Write only
      ReturnMode = ModeWrite;

   } else {
      ABORT_ERROR("IO: Unknown file mode {}: Must be read or write", ModeComp);
   }

   return ReturnMode;

} // End ModeFromString

//------------------------------------------------------------------------------
// Converts string choice for existence behavior to an enum
IfExists IfExistsFromString(
    const std::string &InIfExists // [in] choice of behavior on file existence
) {
   // Set default return value
   IfExists ReturnIfExists = IfExists::Fail;

   // Convert input string to lowercase for easier comparison
   std::string IECompare = InIfExists;
   std::transform(IECompare.begin(), IECompare.end(), IECompare.begin(),
                  [](unsigned char c) { return std::tolower(c); });

   // Determine the appropriate enum to use based on the input string

   if (IECompare == "fail") { // Fail with error
      ReturnIfExists = IfExists::Fail;

   } else if (IECompare == "replace") { // Replace existing file
      ReturnIfExists = IfExists::Replace;

   } else if (IECompare == "append") { // Append or insert into existing file
      ReturnIfExists = IfExists::Append;

   } else {
      ReturnIfExists = IfExists::Fail;
   }

   return ReturnIfExists;

} // End IfExistsFromString

// Methods
//------------------------------------------------------------------------------
// Initializes the IO system based on configuration inputs and
// default MPI communicator
void init(const MPI_Comm &InComm // [in] MPI communicator to use
) {

   // Retrieve parallel IO parameters from the Omega configuration

   Error Err;
   Config *OmegaConfig = Config::getOmegaConfig();

   // Read IO subconfiguration
   Config IOConfig("IO");
   Err = OmegaConfig->get(IOConfig);
   CHECK_ERROR_ABORT(Err, "IO: IO group not found in input Config");

   // Read default file format
   std::string InFileFmt = "netcdf4c"; // set default value
   Err                   = IOConfig.get("IODefaultFormat", InFileFmt);
   CHECK_ERROR_WARN(Err, "IO: DefaultFileFmt not found in Config - using {}",
                    InFileFmt);
   FileFmt DefaultFileFmt = FileFmtFromString(InFileFmt);

   // Read parallel IO settings - default to single-task if config
   // values do not exist
   int NumIOTasks           = 1;
   int IOStride             = 1;
   int IOBaseTask           = 0;
   std::string InRearranger = "box";

   Err = IOConfig.get("IOTasks", NumIOTasks);
   CHECK_ERROR_WARN(Err, "IO: NumIOTasks not found in Config - using {}",
                    NumIOTasks);

   Err = IOConfig.get("IOStride", IOStride);
   CHECK_ERROR_WARN(Err, "IO: IOStride not found in Config - using {}",
                    IOStride);

   Err = IOConfig.get("IOBaseTask", IOBaseTask);
   CHECK_ERROR_WARN(Err, "IO: IOBaseTask not found in Config - using {}",
                    IOBaseTask);

   Err = IOConfig.get("IORearranger", InRearranger);
   CHECK_ERROR_WARN(Err, "IO: Rearranger not found in Config - using {}",
                    InRearranger);
   Rearranger Rearrange = RearrFromString(InRearranger);

   // Call PIO routine to initialize
   DefaultRearr = Rearrange;
   int PIOErr   = PIOc_Init_Intracomm(InComm, NumIOTasks, IOStride, IOBaseTask,
                                      Rearrange, &SysID);
   if (PIOErr != 0)
      ABORT_ERROR("IO::init: Error initializing SCORPIO");

   return;

} // end init

//------------------------------------------------------------------------------
// This routine opens a file for reading or writing, depending on the
// Mode argument. The filename with full path must be supplied and
// a FileID is returned to be used by other IO functions.
// The format of the file is assumed to be the default defined on init
// but can be optionally changed through this open function.
// For files to be written, optional arguments govern the behavior to be
// used if the file already exists, and the precision of any floating point
// variables.
void openFile(
    int &FileID,                 // [out] returned fileID for this file
    const std::string &Filename, // [in] name (incl path) of file to open
    Mode InMode,                 // [in] mode (read or write)
    FileFmt InFormat,            // [in] (optional) file format
    IfExists InIfExists          // [in] (for writes) behavior if file exists
) {

   int PIOErr = 0;        // internal SCORPIO/PIO return call
   int Format = InFormat; // coerce to integer for PIO calls

   switch (InMode) {

   // If reading, open the file for read-only
   case ModeRead:
      PIOErr = PIOc_openfile(SysID, &FileID, &Format, Filename.c_str(), InMode);
      if (PIOErr != PIO_NOERR)
         ABORT_ERROR("IO::openFile: PIO error opening file {} for read",
                     Filename);
      break;

   // If writing, the open will depend on what to do if the file
   // exists.
   case ModeWrite:

      switch (InIfExists) {
      // If the write should be a new file and fail if the
      // file exists, we use create and fail with an error
      case IfExists::Fail:
         PIOErr = PIOc_createfile(SysID, &FileID, &Format, Filename.c_str(),
                                  NC_NOCLOBBER | PIO_64BIT_DATA | InMode);
         if (PIOErr != PIO_NOERR)
            ABORT_ERROR("IO::openFile: PIO error opening file {} for writing",
                        Filename);
         break;

      // If the write should replace any existing file
      // we use create with the CLOBBER option
      case IfExists::Replace:
         PIOErr = PIOc_createfile(SysID, &FileID, &Format, Filename.c_str(),
                                  NC_CLOBBER | PIO_64BIT_DATA | InMode);
         if (PIOErr != PIO_NOERR)
            ABORT_ERROR("IO::openFile: PIO error opening file {} for writing",
                        Filename);
         break;

      // If the write should append or add to an existing file
      // we open the file for reading and writing, but also must create
      // the file if it doesn't already exist
      case IfExists::Append:
         if (std::filesystem::exists(Filename)) {
            PIOErr = PIOc_openfile(SysID, &FileID, &Format, Filename.c_str(),
                                   InMode);
         } else {
            PIOErr = PIOc_createfile(SysID, &FileID, &Format, Filename.c_str(),
                                     PIO_64BIT_DATA | InMode);
         }

         if (PIOErr != PIO_NOERR)
            ABORT_ERROR("IO::openFile: PIO error opening file {} for writing",
                        Filename);
         break;

      default:
         ABORT_ERROR("IO::openFile: unknown IfExists option for writing");

      } // end switch IfExists
      break;

   // Unknown mode
   default:
      ABORT_ERROR("IO::fileOpen: Unknown Mode for file");

   } // End switch on Mode

   return;

} // End openFile
//------------------------------------------------------------------------------
// Closes an open file using the fileID, returns an error code
void closeFile(int &FileID /// [in] ID of the file to be closed
) {
   int Err = 0;

   // Make sure all operations completed before closing
   Err = PIOc_sync(FileID);
   if (Err != PIO_NOERR)
      ABORT_ERROR("IO::closeFile: Error syncing file before closing");

   // Call the PIO close routine
   Err = PIOc_closefile(FileID);
   if (Err != PIO_NOERR)
      ABORT_ERROR("IO::closeFile: Error closing file {} in PIO", FileID);

   return;

} // End closeFile

//------------------------------------------------------------------------------
// Retrieves a dimension length from an input file, given the name
// of the dimension. Returns both the assigned dimension ID and the dimension
// length if exists, but returns an error and a negative length if a dimension
// of that name is not found in the file.
Error getDimFromFile(int FileID, // [in] ID of the file containing dim
                     const std::string &DimName, // [in] name of dimension
                     int &DimID,    // [out] ID assigned to this dimension
                     int &DimLength // [out] global length of the dimension
) {

   Error Err;
   int IOErr = 0;  // Local pio error code
   DimID     = -1; // default (undefined) dimension ID
   DimLength = -1; // default (undefined) global length

   // First get the dimension ID
   IOErr = PIOc_inq_dimid(FileID, DimName.c_str(), &DimID);
   if (IOErr != PIO_NOERR)
      RETURN_ERROR(Err, ErrorCode::Fail,
                   "IO::getDim: PIO error reading dimension {} ", DimName);

   // Now retrieve the length and return
   PIO_Offset InLength;
   IOErr = PIOc_inq_dimlen(FileID, DimID, &InLength);
   if (IOErr == PIO_NOERR) {
      DimLength = InLength;
   } else {
      RETURN_ERROR(Err, ErrorCode::Fail,
                   "IO::getDim: PIO error reading length for dimension {}",
                   DimName);
   }
   return Err;

} // End getDimFromFile

//------------------------------------------------------------------------------
// Defines a dimension for an output file. Returns a dimension id to
// be used in future field definitions as well as an error flag.
int defineDim(int FileID,                 // [in] ID of the file containing dim
              const std::string &DimName, // [in] name of dimension
              int Length                  // [in] length of dimension
) {

   int Err   = 0;
   int DimID = -1;

   Err = PIOc_def_dim(FileID, DimName.c_str(), Length, &DimID);
   if (Err != PIO_NOERR) {
      ABORT_ERROR("IO::defineDim: PIO error while defining dimension {}",
                  DimName);
   }

   return DimID;

} // End defineDim

//------------------------------------------------------------------------------
// Writes metadata (name, value) associated with a variable and/or file.
// The variable ID can be GlobalID for global file and/or simulation
// metadata. Specific interfaces for each supported data type are mapped to
// the generic interface rather than templating, providing a cleaner interface.
void writeMeta(const std::string &MetaName, // [in] name of metadata
               I4 MetaValue,                // [in] value of metadata
               int FileID,                  // [in] ID of the file for writing
               int VarID // [in] ID for variable associated with metadata
) {

   int Err             = 0;
   IODataType MetaType = IOTypeI4;
   PIO_Offset Length   = 1;

   // Check to see if metadata already exists and has the same value
   I4 TmpValue;
   Err = PIOc_get_att(FileID, VarID, MetaName.c_str(), &TmpValue);
   if (Err == PIO_NOERR) { // Metadata already exists, check value is same
      if (TmpValue == MetaValue)
         return;
   }

   // Write the metadata
   Err = PIOc_put_att(FileID, VarID, MetaName.c_str(), MetaType, Length,
                      &MetaValue);
   if (Err != PIO_NOERR)
      ABORT_ERROR("IO::writeMeta(I4): PIO error while writing metadata {}",
                  MetaName);

   return;

} // End writeMeta (I4)

void writeMeta(const std::string &MetaName, // [in] name of metadata
               I8 MetaValue,                // [in] value of metadata
               int FileID,                  // [in] ID of the file for writing
               int VarID // [in] ID for variable associated with metadata
) {

   int Err             = 0;
   IODataType MetaType = IOTypeI8;
   PIO_Offset Length   = 1;

   // Check to see if metadata already exists and has the same value
   I8 TmpValue;
   Err = PIOc_get_att(FileID, VarID, MetaName.c_str(), &TmpValue);
   if (Err == PIO_NOERR) { // Metadata already exists, if value is same return
      if (TmpValue == MetaValue)
         return;
   }

   // Write the metadata
   Err = PIOc_put_att(FileID, VarID, MetaName.c_str(), MetaType, Length,
                      &MetaValue);
   if (Err != PIO_NOERR)
      ABORT_ERROR("IO::writeMeta(I8): PIO error while writing metadata {}",
                  MetaName);

   return;

} // End writeMeta (I8)

void writeMeta(const std::string &MetaName, // [in] name of metadata
               R4 MetaValue,                // [in] value of metadata
               int FileID,                  // [in] ID of the file for writing
               int VarID // [in] ID for variable associated with metadata
) {

   int Err             = 0;
   IODataType MetaType = IOTypeR4;
   PIO_Offset Length   = 1;

   // Check to see if metadata already exists and has the same value
   R4 TmpValue;
   Err = PIOc_get_att(FileID, VarID, MetaName.c_str(), &TmpValue);
   if (Err == PIO_NOERR) { // Metadata already exists, check value is same
      if (TmpValue == MetaValue)
         return;
   }

   // Write the metadata
   Err = PIOc_put_att(FileID, VarID, MetaName.c_str(), MetaType, Length,
                      &MetaValue);
   if (Err != PIO_NOERR)
      ABORT_ERROR("IO::writeMeta(R4): PIO error while writing metadata {}",
                  MetaName);

   return;

} // End writeMeta (R4)

void writeMeta(const std::string &MetaName, // [in] name of metadata
               R8 MetaValue,                // [in] value of metadata
               int FileID,                  // [in] ID of the file for writing
               int VarID // [in] ID for variable associated with metadata
) {

   int Err             = 0;
   IODataType MetaType = IOTypeR8;
   PIO_Offset Length   = 1;

   // Check to see if metadata already exists and has the same value
   R8 TmpValue;
   Err = PIOc_get_att(FileID, VarID, MetaName.c_str(), &TmpValue);
   if (Err == PIO_NOERR) { // Metadata already exists, check value is same
      if (TmpValue == MetaValue)
         return;
   }

   // Write the metadata
   Err = PIOc_put_att(FileID, VarID, MetaName.c_str(), MetaType, Length,
                      &MetaValue);
   if (Err != PIO_NOERR)
      ABORT_ERROR("IO::writeMeta(R8): PIO error while writing metadata {}",
                  MetaName);

   return;

} // End writeMeta (R8)

void writeMeta(const std::string &MetaName,  // [in] name of metadata
               const std::string &MetaValue, // [in] value of metadata
               int FileID,                   // [in] ID of the file for writing
               int VarID // [in] ID for variable associated with metadata
) {

   int Err             = 0;
   IODataType MetaType = IOTypeChar;
   PIO_Offset Length   = MetaValue.length() + 1; // add 1 for char terminator

   // Check to see if metadata already exists and has the same value
   // For strings, we need to get the length first to make sure the
   // string is long enough to read the value
   PIO_Offset TmpLength;
   Err = PIOc_inq_attlen(FileID, VarID, MetaName.c_str(), &TmpLength);
   if (Err == PIO_NOERR and TmpLength > 0) { // Metadata exists and has length
      std::string TmpValue;
      TmpValue.resize(TmpLength);
      Err = PIOc_get_att(FileID, VarID, MetaName.c_str(),
                         (void *)TmpValue.c_str());
      if (Err == PIO_NOERR) {
         if (TmpValue == MetaValue)
            return;
      }
   }

   // Write the metadata
   Err = PIOc_put_att(FileID, VarID, MetaName.c_str(), MetaType, Length,
                      (void *)MetaValue.c_str());
   if (Err != PIO_NOERR)
      ABORT_ERROR("IO::writeMeta(str): PIO error while writing metadata {}",
                  MetaName);

   return;

} // End writeMeta (str)

//------------------------------------------------------------------------------
// Reads metadata (name, value) associated with a variable and/or file.
// The variable ID can be GlobalID for global file and/or simulation
// metadata. Specific interfaces for supported data types are aliased to the
// same generic form. The functions also return an error code so the
// calling routine can decide actions if the metadata is not found.
Error readMeta(const std::string &MetaName, // [in] name of metadata
               I4 &MetaValue,               // [out] value of metadata
               int FileID,                  // [in] ID of the file for writing
               int VarID // [in] ID for variable associated with metadata
) {
   Error Err;
   int PIOErr = 0;

   PIOErr = PIOc_get_att(FileID, VarID, MetaName.c_str(), &MetaValue);
   if (PIOErr != PIO_NOERR)
      RETURN_ERROR(Err, ErrorCode::Fail,
                   "IO::readMeta(I4): PIO error while reading metadata for {}",
                   MetaName);

   return Err;

} // End readMeta (I4)

Error readMeta(const std::string &MetaName, // [in] name of metadata
               I8 &MetaValue,               // [out] value of metadata
               int FileID,                  // [in] ID of the file for writing
               int VarID // [in] ID for variable associated with metadata
) {
   Error Err;
   int PIOErr = 0;

   PIOErr = PIOc_get_att(FileID, VarID, MetaName.c_str(), &MetaValue);
   if (PIOErr != PIO_NOERR)
      RETURN_ERROR(Err, ErrorCode::Fail,
                   "IO::readMeta(I8): PIO error while reading metadata for {}",
                   MetaName);

   return Err;

} // End readMeta (I8)

Error readMeta(const std::string &MetaName, // [in] name of metadata
               R4 &MetaValue,               // [out] value of metadata
               int FileID,                  // [in] ID of the file for writing
               int VarID // [in] ID for variable associated with metadata
) {
   Error Err;
   int PIOErr = 0;

   PIOErr = PIOc_get_att(FileID, VarID, MetaName.c_str(), &MetaValue);
   if (PIOErr != PIO_NOERR)
      RETURN_ERROR(Err, ErrorCode::Fail,
                   "IO::readMeta(R4): PIO error while reading metadata for {}",
                   MetaName);

   return Err;

} // End readMeta (R4)

Error readMeta(const std::string &MetaName, // [in] name of metadata
               R8 &MetaValue,               // [out] value of metadata
               int FileID,                  // [in] ID of the file for writing
               int VarID // [in] ID for variable associated with metadata
) {
   Error Err;
   int PIOErr = 0;

   PIOErr = PIOc_get_att(FileID, VarID, MetaName.c_str(), &MetaValue);
   if (PIOErr != PIO_NOERR)
      RETURN_ERROR(Err, ErrorCode::Fail,
                   "IO::readMeta(R8): PIO error while reading metadata for {}",
                   MetaName);

   return Err;

} // End readMeta (R8)

Error readMeta(const std::string &MetaName, // [in] name of metadata
               std::string &MetaValue,      // [out] value of metadata
               int FileID,                  // [in] ID of the file for writing
               int VarID // [in] ID for variable associated with metadata
) {
   Error Err;
   int PIOErr = 0;

   // For string variables, find the length of the string first
   // and resize the string to make sure the allocated length is large
   // enough when passing the pointer later
   PIO_Offset Length;
   PIOErr = PIOc_inq_attlen(FileID, VarID, MetaName.c_str(), &Length);
   if (PIOErr != PIO_NOERR)
      RETURN_ERROR(Err, ErrorCode::Fail,
                   "IO::readMeta(str): PIO error while finding string length");

   int StrLength = Length - 1;
   MetaValue.resize(StrLength);

   // Now read the string
   PIOErr =
       PIOc_get_att(FileID, VarID, MetaName.c_str(), (void *)MetaValue.c_str());
   if (PIOErr != PIO_NOERR)
      RETURN_ERROR(Err, ErrorCode::Fail,
                   "IO::readMeta(str): PIO error while reading metadata for {}",
                   MetaName);

   return Err;

} // End readMeta (str)

//------------------------------------------------------------------------------
// Defines a variable for an output file. The name and dimensions of
// the variable must be supplied. An ID is assigned to the variable and returned
// for later use in the writing of the variable.
int defineVar(int FileID,                 // [in] ID of the file containing var
              const std::string &VarName, // [in] name of variable
              IODataType VarType,         // [in] data type for the variable
              int NDims,                  // [in] number of dimensions
              int *DimIDs                 // [in] vector of NDims dimension IDs
) {

   int PIOErr = 0;
   int VarID  = -1;

   // First check to see if the variable exists (if reading or if appending
   // to an existing file)
   PIOErr = PIOc_inq_varid(FileID, VarName.c_str(), &VarID);

   // If the variable is not found, define the new variable
   if (PIOErr != PIO_NOERR) {

      PIOErr =
          PIOc_def_var(FileID, VarName.c_str(), VarType, NDims, DimIDs, &VarID);
      if (PIOErr != PIO_NOERR)
         ABORT_ERROR("IO::defineVar: PIO error while defining variable {}",
                     VarName);
   }

   return VarID;

} // End defineVar

//------------------------------------------------------------------------------
/// Ends define phase signifying all field definitions and metadata
/// have been written and the larger data sets can now be written
void endDefinePhase(int FileID ///< [in] ID of the file being written
) {
   int Err = 0;

   Err = PIOc_enddef(FileID);
   if (Err != PIO_NOERR)
      ABORT_ERROR("IO::endDefinePhase: PIO error in enddef");

} // End endDefinePhase

//------------------------------------------------------------------------------
// Creates a PIO decomposition description to describe the layout of
// a distributed array of given type. It returns an ID for future use.
int createDecomp(
    IODataType VarType,                 // [in] data type of array
    int NDims,                          // [in] number of array dimensions
    const std::vector<int> &DimLengths, // [in] global dimension lengths
    int Size,                           // [in] local size of array
    const std::vector<int> &GlobalIndx, // [in] global indx for each local indx
    Rearranger Rearr                    // [in] rearranger method to use
) {

   int DecompID = -1;
   int PIOErr   = 0; // internal error code for PIO calls

   // Convert global index array into an offset array expected by PIO
   std::vector<PIO_Offset> CompMap;
   CompMap.resize(Size);
   for (int i = 0; i < Size; ++i) {
      CompMap[i] = GlobalIndx[i];
   } // end global index loop

   // Call the PIO routine to define the decomposition
   // int TmpRearr = Rearr; // needed for type compliance across interface
   // Err = PIOc_InitDecomp(SysID, VarType, NDims, DimLengths, Size, CompMap,
   //                      &DecompID, &TmpRearr, nullptr, nullptr);
   PIOErr = PIOc_init_decomp(SysID, VarType, NDims, &DimLengths[0], Size,
                             &CompMap[0], &DecompID, Rearr, nullptr, nullptr);
   if (PIOErr != PIO_NOERR)
      ABORT_ERROR("IO::createDecomp: PIO error defining decomposition");

   return DecompID;

} // End createDecomp

//------------------------------------------------------------------------------
// Removes a defined PIO decomposition description to free memory
void destroyDecomp(int &DecompID // [inout] ID for decomposition to be removed
) {

   int Err = PIOc_freedecomp(SysID, DecompID);
   if (Err != PIO_NOERR)
      ABORT_ERROR("IO::destroyDecomp: PIO error freeing decomposition");

} // End destroyDecomp

//------------------------------------------------------------------------------
// Reads a distributed array. Uses a void pointer for generic interface.
// All arrays are assumed to be in contiguous storage. Returns an error code
// that is mostly for use by the calling routine to retry reads.
Error readArray(void *Array,                // [out] array to be read
                int Size,                   // [in] local size of array
                const std::string &VarName, // [in] name of variable to read
                int FileID,                 // [in] ID of open file to read from
                int DecompID,               // [in] decomposition ID for var
                int &VarID, // [out] Id assigned to variable for later use
                int Frame   // [in] opt frame if multiple time slices
) {

   Error Err;      // returned error code
   int PIOErr = 0; // internal error code for PIO calls

   // Find variable ID from file
   PIOErr = PIOc_inq_varid(FileID, VarName.c_str(), &VarID);
   if (PIOErr != PIO_NOERR)
      RETURN_ERROR(Err, ErrorCode::Fail,
                   "IO::readArray: Error finding varid for variable {}",
                   VarName);

   if (Frame >= 0) {
      PIOErr = PIOc_setframe(FileID, VarID, Frame);
      if (PIOErr != PIO_NOERR)
         RETURN_ERROR(Err, ErrorCode::Fail,
                      "Error setting frame while reading distributed array");
   }

   // PIO Read array call to read the distributed array
   PIO_Offset ASize = Size;
   PIOErr           = PIOc_read_darray(FileID, VarID, DecompID, ASize, Array);
   if (PIOErr != PIO_NOERR)
      RETURN_ERROR(Err, ErrorCode::Fail,
                   "IO::readArray: Error in SCORPIO read array for variable {}",
                   VarName);

   return Err;

} // End IOReadArray

//------------------------------------------------------------------------------
// Reads a non-distributed variable. Uses a void pointer for generic interface.
// All arrays are assumed to be in contiguous storage. Returns an error code so
// that the calling routine can re-try on failure.
Error readNDVar(void *Variable,             // [out] array to be read
                const std::string &VarName, // [in] name of variable to read
                int FileID,                 // [in] ID of open file to read from
                int &VarID, // [out] Id assigned to variable for later use
                int Frame,  // [in] opt frame/slice for multiframe streams
                std::vector<int> *DimLengths // [in] dim lengths for multiframe
) {

   Error Err;      // returned error code
   int PIOErr = 0; // internal error code for PIO calls

   // Find variable ID from file
   PIOErr = PIOc_inq_varid(FileID, VarName.c_str(), &VarID);
   if (PIOErr != PIO_NOERR)
      RETURN_ERROR(Err, ErrorCode::Fail,
                   "IO::readArray: Error finding varid for variable {}",
                   VarName);

   if (Frame >= 0) { // time dependent field so must use get_vara

      int NDims = DimLengths->size();
      // Start and count arguments include the unlimited time dim (Frame)
      std::vector<PIO_Offset> Start(NDims + 1, 0);
      std::vector<PIO_Offset> Count(NDims + 1, 0);
      Start[0] = Frame;
      Count[0] = 1;
      for (int IDim = 1; IDim <= NDims; ++IDim) {
         Count[IDim] = DimLengths->at(IDim - 1);
      }

      PIOErr =
          PIOc_get_vara(FileID, VarID, Start.data(), Count.data(), Variable);
      if (PIOErr != PIO_NOERR)
         RETURN_ERROR(Err, ErrorCode::Fail,
                      "IO::readNDVar: Error in PIO get_vara for variable {}",
                      VarName);

   } else { // Not a time-dependent field, can use default get_var

      // PIO get call to read non-distributed array
      PIOErr = PIOc_get_var(FileID, VarID, Variable);
      if (PIOErr != PIO_NOERR)
         RETURN_ERROR(Err, ErrorCode::Fail,
                      "IO::readNDVar: Error in PIO get_var for variable {}",
                      VarName);
   }

   return Err;

} // End IOReadNDVar

//------------------------------------------------------------------------------
// Writes a distributed array. This generic interface uses void pointers.
// All arrays are assumed to be in contiguous storage and the variable
// must have a valid ID assigned by the defineVar function.

void writeArray(void *Array,     // [in] array to be written
                int Size,        // [in] size of array to be written
                void *FillValue, // [in] value to use for missing entries
                int FileID,      // [in] ID of open file to write to
                int DecompID,    // [in] decomposition ID for this var
                int VarID,       // [in] variable ID assigned by defineVar
                int Frame        // [in] frame/slice for multiframe streams
) {
   int Err = 0; // internal error code for PIO calls

   if (Frame >= 0) {
      Err = PIOc_setframe(FileID, VarID, Frame);
      if (Err != PIO_NOERR)
         ABORT_ERROR("IO::writeArray: Error setting frame");
   }

   PIO_Offset Asize = Size;

   Err = PIOc_write_darray(FileID, VarID, DecompID, Asize, Array, FillValue);
   if (Err != PIO_NOERR)
      ABORT_ERROR("IO::writeArray: Error in PIO writing distributed array");

   // Make sure write is complete before returning
   // We may be able to remove this for efficiency later but it was
   // needed during testing
   Err = PIOc_sync(FileID);
   if (Err != PIO_NOERR)
      LOG_WARN("IO::writeArray: Error in PIO sychronizing file after write");

   return;

} // end writeArray

//------------------------------------------------------------------------------
// Writes a non-distributed variable. This generic interface uses void pointers.
// All arrays are assumed to be in contiguous storage and the variable
// must have a valid ID assigned by the defineVar function. For time-dependent
// fields, optional arguments for the frame of the unlimited time axis and a
// vector of dimension lengths must be provided.

void writeNDVar(void *Variable, // [in] variable to be written
                int FileID,     // [in] ID of open file to write to
                int VarID,      // [in] variable ID assigned by defineVar
                int Frame,      // [in] frame/slice for multiframe streams
                std::vector<int> *DimLengths // [in] dimension lengths
) {
   int Err = 0; // internal error code for PIO calls

   if (Frame >= 0) { // time dependent field so must use put_vara
      int NDims = DimLengths->size();
      // Start and count arguments include the unlimited time dim (Frame)
      std::vector<PIO_Offset> Start(NDims + 1, 0);
      std::vector<PIO_Offset> Count(NDims + 1, 0);
      Start[0] = Frame;
      Count[0] = 1;
      for (int IDim = 1; IDim <= NDims; ++IDim) {
         Count[IDim] = DimLengths->at(IDim - 1);
      }

      Err = PIOc_put_vara(FileID, VarID, Start.data(), Count.data(), Variable);
      if (Err != PIO_NOERR)
         ABORT_ERROR("IO::writeNDVar: Error in PIO put_vara while writing");

   } else {

      Err = PIOc_put_var(FileID, VarID, Variable);
      if (Err != PIO_NOERR)
         ABORT_ERROR("IO::writeNDVar: Error in PIO put_var while writing");
   }

   return;

} // end writeNDVar

//------------------------------------------------------------------------------

} // end namespace IO
} // end namespace OMEGA

//===----------------------------------------------------------------------===//
