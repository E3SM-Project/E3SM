//===-- infra/IOStream.cpp - IO stream implementation -----------*- C++ -*-===//
//
// This file implements classes and methods for IO Streams. IOStreams define
// reading and writing of fields from/to a data file. Each stream contains
// information on the file location, the frequency of input/output, the
// fields to be read/written and other information necessary for file I/O.
// Note that this is different from the C++ stdlib iostreams
//
//===----------------------------------------------------------------------===//

#include "IOStream.h"
#include "DataTypes.h"
#include "Dimension.h"
#include "Error.h"
#include "Field.h"
#include "IO.h"
#include "Logging.h"
#include "OmegaKokkos.h"
#include "TimeMgr.h"
#include <algorithm>
#include <any>
#include <cctype>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <typeindex>
#include <typeinfo>

namespace OMEGA {

// Create static class members
std::map<std::string, std::shared_ptr<IOStream>> IOStream::AllStreams;

//------------------------------------------------------------------------------
// Initializes all streams defined in the input configuration file. This
// does not validate the contents of the streams since the relevant Fields
// may not have been defined yet.
void IOStream::init(Clock *&ModelClock //< [inout] Omega model clock
) {

   Error Err;

   // Retrieve the model configuration and get the streams sub-config
   Config *OmegaConfig = Config::getOmegaConfig();
   Config StreamsCfgAll("IOStreams");
   Err = OmegaConfig->get(StreamsCfgAll);
   CHECK_ERROR_ABORT(
       Err, "Could not find Streams configuration in Omega input file");

   // Loop over all streams in the subconfiguration and create them
   for (auto It = StreamsCfgAll.begin(); It != StreamsCfgAll.end(); ++It) {

      // Get the stream name and the sub-configuration associated with it
      std::string StreamName = It->first.as<std::string>();
      Config StreamCfg(StreamName);
      Err = StreamsCfgAll.get(StreamCfg);
      CHECK_ERROR_ABORT(Err, "Error retrieving configuration for stream {}",
                        StreamName);

      // Call the create routine to create the stream
      create(StreamName, StreamCfg, ModelClock);

   } // end loop over all streams

} // End initialize

//------------------------------------------------------------------------------
/// Initialize with no clock. Convenience overload for testing.
void IOStream::init(void) {

   // Create an empty dummy clock
   Clock *ModelClock = new Clock;
   init(ModelClock);

} // End initialize(void)

//------------------------------------------------------------------------------
// Performs a final write of any streams that have the OnShutdown option and
// then removes all streams to clean up.
void IOStream::finalize(
    const Clock *ModelClock // [in] Model clock needed for time stamps
) {

   bool FinalCall = true;

   // Loop over all streams and call write function for any write streams
   // with the OnShutdown flag
   for (auto Iter = AllStreams.begin(); Iter != AllStreams.end(); Iter++) {

      std::string StreamName               = Iter->first;
      std::shared_ptr<IOStream> ThisStream = Iter->second;
      bool ForceWrite                      = false;

      if (ThisStream->OnShutdown)
         ThisStream->writeStream(ModelClock, ForceWrite, FinalCall);

   } // end loop over streams

   // Remove all streams
   AllStreams.clear();

   return;

} // End finalize

//------------------------------------------------------------------------------
// Finalize with no clock. Convenience overload for testing.
void IOStream::finalize(void) {

   // Create an empty dummy clock
   Clock *ModelClock = new Clock;
   finalize(ModelClock);

} // End finalize(void)

//------------------------------------------------------------------------------
// Retrieves a pointer to a previously defined stream.
std::shared_ptr<IOStream>
IOStream::get(const std::string &StreamName ///< [in] name of stream to retrieve
) {
   // Check that stream exists
   if (AllStreams.find(StreamName) == AllStreams.end())
      ABORT_ERROR("IOStream: Cannot retrieve Stream {}. Stream does not exist.",
                  StreamName);

   // Return the pointer from the map of all streams
   return AllStreams[StreamName];

} // End get stream

//------------------------------------------------------------------------------
// Adds a field to the contents of a stream. Because streams may be created
// before all Fields have been defined, we only store the name. Validity
// is either checked during read/write or can be checked using the validate
// function.
void IOStream::addField(const std::string &FieldName ///< [in] Name of field
) {
   this->Contents.insert(FieldName);
} // End addField

//------------------------------------------------------------------------------
// Removes a field from the contents. Provided for symmetry, but not
// typically used.
void IOStream::removeField(const std::string &FieldName ///< [in] Name of field
) {
   this->Contents.erase(FieldName);
} // End removeField

//------------------------------------------------------------------------------
// Validates an IOStream and its contents. This must be called after all
// relevant fields have been defined, typically toward the end of omega model
// initialization. The routine also expands all group names to the individual
// field names that are members of the group. Returns true if all contents and
// variables are valid.
bool IOStream::validate() {

   bool ReturnVal = true;
   // Return if already validated
   if (Validated)
      return ReturnVal;

   // Expand group names to list of individual fields
   // First identify any group names in the Contents
   std::set<std::string> GroupNames;
   for (auto IField = Contents.begin(); IField != Contents.end(); ++IField) {
      std::string FieldName = *IField;
      // If this name is a group name, add it to the list
      if (FieldGroup::exists(FieldName))
         GroupNames.insert(FieldName);
   }

   // Now for each group, extract field names and add to contents if it
   // is not there already. Since Contents is a set we can just insert
   // and it will take care of duplicate fields.
   for (auto IGrp = GroupNames.begin(); IGrp != GroupNames.end(); ++IGrp) {
      std::string GrpName = *IGrp;
      // Get Fields from Group
      std::set<std::string> FieldList =
          FieldGroup::getFieldListFromGroup(GrpName);
      // Insert field names into Contents
      for (auto IField = FieldList.begin(); IField != FieldList.end();
           ++IField) {
         Contents.insert(*IField);
      }
      // Remove group name from contents
      Contents.erase(GrpName);
   }

   // Loop through all the field names in Contents and check whether they
   // have been defined as a Field
   for (auto IField = Contents.begin(); IField != Contents.end(); ++IField) {
      std::string FieldName = *IField;

      if (not Field::exists(FieldName))
         ABORT_ERROR("IOStream validate: "
                     "Cannot validate stream {}: Field {} has not been defined",
                     Name, FieldName);
   }

   if (ReturnVal)
      Validated = true;
   return ReturnVal;

} // End validate

//------------------------------------------------------------------------------
// Checks that a stream has been validated
bool IOStream::isValidated() { return Validated; }

//------------------------------------------------------------------------------
// Validates all streams and their contents. If used, this must be called at
// the end of initialization to ensure all Fields have been defined.
// Returns true if all streams are valid.
bool IOStream::validateAll() {

   bool ReturnVal = true; // default is all valid

   // Loop over all streams and call validate function
   for (auto Iter = AllStreams.begin(); Iter != AllStreams.end(); Iter++) {
      std::string StreamName     = Iter->first;
      std::shared_ptr ThisStream = Iter->second;
      bool Valid                 = ThisStream->validate();
      if (!Valid) {
         ReturnVal = false;
         ABORT_ERROR("IOStream validateAll: stream {} has invalid entries",
                     StreamName);
      }
   }

   return ReturnVal;

} // End validateAll

//------------------------------------------------------------------------------
// Reads a single stream if it is time.
Error IOStream::read(
    const std::string &StreamName, // [in] Name of stream
    const Clock *ModelClock,       // [in] Model clock for time info
    Metadata &ReqMetadata, // [inout] global metadata requested from file
    bool ForceRead         // [in] optional: read even if not time
) {

   Error Err; // returned error code

   // Retrieve stream by name and make sure it exists
   auto StreamItr = AllStreams.find(StreamName);
   if (StreamItr == AllStreams.end())
      RETURN_ERROR(Err, ErrorCode::Fail, "IOStream::Read: Stream {} not found",
                   StreamName);

   // Stream found, call the read function
   std::shared_ptr<IOStream> ThisStream = StreamItr->second;
   Err += ThisStream->readStream(ModelClock, ReqMetadata, ForceRead);

   return Err;

} // End read stream

//------------------------------------------------------------------------------
// Reads a single stream regardless of time. Returns an error code.
Error IOStream::read(const std::string &StreamName // [in] Name of stream
) {

   Error Err; // returned error code

   // create empty Clock and Metadata
   Clock *ModelClock = new Clock;
   Metadata ReqMetaData;
   bool ForceRead = true;

   Err = read(StreamName, ModelClock, ReqMetaData, ForceRead);

   return Err;

} // End read(void)

//------------------------------------------------------------------------------
// Writes a single stream if it is time.
void IOStream::write(
    const std::string &StreamName, // [in] Name of stream
    const Clock *ModelClock,       // [in] Model clock needed for time stamps
    bool ForceWrite                // [in] optional: write even if not time
) {

   // Retrieve stream by name and make sure it exists
   auto StreamItr = AllStreams.find(StreamName);
   if (StreamItr == AllStreams.end())
      ABORT_ERROR("Unable to write stream {}. Stream not defined", StreamName);

   // Stream found, call the write function
   std::shared_ptr<IOStream> ThisStream = StreamItr->second;
   ThisStream->writeStream(ModelClock, ForceWrite);

   return;

} // End write stream

//------------------------------------------------------------------------------
// Loops through all streams and writes them if it is time. This is
// useful if most I/O is consolidated at one point (eg end of step).
void IOStream::writeAll(
    const Clock *ModelClock // [in] Model clock needed for time stamps
) {

   // Loop over all streams and call write function for any write streams
   for (auto Iter = AllStreams.begin(); Iter != AllStreams.end(); Iter++) {

      std::string StreamName               = Iter->first;
      std::shared_ptr<IOStream> ThisStream = Iter->second;

      if (ThisStream->Mode == IO::ModeWrite) {
         ThisStream->writeStream(ModelClock);
      }
   }

   return;

} // End writeAll

//------------------------------------------------------------------------------
// Constructs an empty IOStream
IOStream::IOStream() {

   Name               = "Unknown";
   Filename           = "Unknown";
   FilenameIsTemplate = false;
   ExistAction        = IO::IfExists::Fail;
   Mode               = IO::Mode::ModeUnknown;
   ReducePrecision    = false;
   OnStartup          = false;
   OnShutdown         = false;
   UsePointer         = false;
   PtrFilename        = " ";
   UseStartEnd        = false;
   Validated          = false;
}

//------------------------------------------------------------------------------
// Creates a new stream and adds to the list of all streams, based on
// options in the input model configuration. This routine is called by
// the IOStreams initialize function. It requires an initialized model
// clock and stream alarms are attached to this clock during creation.

void IOStream::create(const std::string &StreamName, //< [in] name of stream
                      Config &StreamConfig, //< [in] input stream configuration
                      Clock *&ModelClock    //< [inout] Omega model clock
) {

   Error Err; // internal error code

   // Check whether the stream already exists
   if (AllStreams.find(StreamName) != AllStreams.end())
      ABORT_ERROR("IOStream: Attempt to create stream {} that already exists",
                  StreamName);

   // Create a new pointer and set name
   auto NewStream  = std::make_shared<IOStream>();
   NewStream->Name = StreamName;

   // Set file mode (Read/Write)
   std::string StreamMode;
   Err += StreamConfig.get("Mode", StreamMode);
   NewStream->Mode = IO::ModeFromString(StreamMode);
   if (Err.isFail() or (NewStream->Mode == IO::ModeUnknown))
      ABORT_ERROR("IOStream::create: Bad or non-existent Mode for IO stream {}",
                  StreamName);

   // Set file name
   // First check if using a pointer file
   Err += StreamConfig.get("UsePointerFile", NewStream->UsePointer);
   CHECK_ERROR_ABORT(Err, "UsePointerFile flag not found in IO stream {}",
                     StreamName);

   // If the stream uses a pointer file, get the pointer filename
   if (NewStream->UsePointer) {
      NewStream->PtrFilename = ""; // initialize to blank filename
      Err += StreamConfig.get("PointerFilename", NewStream->PtrFilename);
      CHECK_ERROR_ABORT(Err, "Pointer filename not found for stream {}",
                        StreamName);
   }

   // For file reads that are using a pointer file, the filename is
   // read later from that pointer file. All other cases need to read
   // a filename or filename template from config
   // Otherwise, read filename from config
   std::string StreamFilename = "";
   if (NewStream->Mode != IO::ModeRead or !(NewStream->UsePointer)) {
      Err += StreamConfig.get("Filename", StreamFilename);
      CHECK_ERROR_ABORT(Err,
                        "Error getting Filename for IO stream {} from Config",
                        StreamName);

      // Check to see if filename is a template and needs to be constructed
      // later with time information
      auto TemplateFound = StreamFilename.find("$");
      if (TemplateFound != std::string::npos) {
         NewStream->FilenameIsTemplate = true;
      } else {
         NewStream->FilenameIsTemplate = false;
      }
   }
   // Add filename to stream
   NewStream->Filename = StreamFilename;

   // Set flag to reduce precision for double precision reals. If no flag
   // present, assume full (double) precision
   std::string PrecisionString;
   Err += StreamConfig.get("Precision", PrecisionString);
   if (Err.isFail())
      PrecisionString = "double";
   NewStream->setPrecisionFlag(PrecisionString);

   // Set the action to take if a file already exists
   // This is only needed for writes so only perform check for write mode
   NewStream->ExistAction = IO::IfExists::Fail; // default is to fail
   if (NewStream->Mode == IO::ModeWrite) {
      std::string ExistAct;
      Err += StreamConfig.get("IfExists", ExistAct);
      if (Err.isSuccess())
         NewStream->ExistAction = IO::IfExistsFromString(ExistAct);
   }

   // Set alarm based on read/write frequency
   // Use stream name as alarm name
   std::string AlarmName = StreamName;

   // For alarms, need to retrieve clock start time
   TimeInstant ClockStart = ModelClock->getStartTime();

   // Read frequency of input/output
   int IOFreq;
   Err += StreamConfig.get("Freq", IOFreq);
   CHECK_ERROR_ABORT(Err, "Frequency missing for stream {}", StreamName);

   if (IOFreq < 1) {
      ABORT_ERROR("Invalid frequency {} for IO stream {} ", StreamName);
   }
   std::string IOFreqUnits;
   Err += StreamConfig.get("FreqUnits", IOFreqUnits);
   CHECK_ERROR_ABORT(Err, "FreqUnits missing for stream {}", StreamName);

   // convert string to lower case for easier comparison
   std::transform(IOFreqUnits.begin(), IOFreqUnits.end(), IOFreqUnits.begin(),
                  [](unsigned char C) { return std::tolower(C); });

   // Based on input frequency and units, create the alarm or set flags
   bool HasAlarm         = false;
   NewStream->OnStartup  = false;
   NewStream->OnShutdown = false;

   // Convert to standard freq units if possible
   TimeUnits IOFrqUnits = TimeUnitsFromString(IOFreqUnits);
   if (IOFrqUnits != TimeUnits::None) { // standard time units cases

      TimeInterval AlarmInt(IOFreq, IOFrqUnits);
      NewStream->MyAlarm = Alarm(AlarmName, AlarmInt, ClockStart);
      HasAlarm           = true;

   } else if (IOFreqUnits == "onstartup") { // special startup case

      NewStream->OnStartup = true;

   } else if (IOFreqUnits == "onshutdown") { // special shutdown case

      NewStream->OnShutdown = true;

      // single time instant case
   } else if (IOFreqUnits == "attime" or IOFreqUnits == "ontime" or
              IOFreqUnits == "time" or IOFreqUnits == "timeinstant") {

      // A one-time event for this stream - use the StartTime string
      // as the time instant to use
      std::string StrtTime;
      Err += StreamConfig.get("StartTime", StrtTime);
      if (Err.isSuccess()) {
         TimeInstant AlarmTime(StrtTime);
         NewStream->MyAlarm = Alarm(AlarmName, AlarmTime);
         HasAlarm           = true;
      } else {
         ABORT_ERROR("Stream {} requests a one-time read/write but StartTime"
                     "not provided",
                     StreamName);
      }

   } else if (IOFreqUnits == "never") { // special never case

      LOG_WARN("Stream {} has IO frequency of never and will be skipped",
               StreamName);
      return;

   } else { // default error

      if (!NewStream->OnStartup and !NewStream->OnShutdown) {
         ABORT_ERROR("Unknown IOFreqUnits option for stream {}", StreamName);
      }

   } // end if IOFreqUnits
   // If an alarm is set, attach it to the model clock
   if (HasAlarm)
      ModelClock->attachAlarm(&(NewStream->MyAlarm));

   // Check the file frequency to see if there will be multiple records or
   // time slices per file
   NewStream->Multiframe = false;
   int FileFreq          = 0;
   std::string FileUnitStr;

   Error Err1          = StreamConfig.get("FileFreq", FileFreq);
   Error Err2          = StreamConfig.get("FileFreqUnits", FileUnitStr);
   TimeUnits FileUnits = TimeUnitsFromString(FileUnitStr);

   // Various error checks
   if (Err1.isSuccess() or Err2.isSuccess()) {
      // file freq variable present, assume Multiframe and check validity
      NewStream->Multiframe = true;

      if (!Err1.isSuccess()) {
         ABORT_ERROR("File frequency units provided but file freq missing "
                     "for stream {}",
                     StreamName);
      }
      if (!Err2.isSuccess()) {
         ABORT_ERROR("File frequency provided but file freq units missing "
                     "for stream {}",
                     StreamName);
      }
      if (FileFreq < 1) {
         ABORT_ERROR("Invalid file frequency {} for IO stream {} ", StreamName);
      }
      if (FileUnits == TimeUnits::None) {
         ABORT_ERROR("Bad units for file frequency in IO stream {} ",
                     StreamName);
      }
      if (NewStream->Mode == IO::ModeRead) {
         ABORT_ERROR("Multiple time slices currently not supported for input "
                     "stream {} ",
                     StreamName);
      }
   }

   // Create and attach alarm for file frequency if needed
   if (NewStream->Multiframe) {

      // Create time interval for files and check compatibility
      // with data IO frequency
      TimeInterval FileInt(FileFreq, FileUnits);

      // Cannot compare calendar and non-calendar intervals, so need
      // need a better check here for incompatible intervals
      // (eg file freq greater than data freq)

      // create and attach alarm
      std::string FileAlarmName = StreamName + "File";
      NewStream->FileAlarm      = Alarm(FileAlarmName, FileInt, ClockStart);
      ModelClock->attachAlarm(&(NewStream->FileAlarm));

   } // end multiframe

   // Use a start and end time to define an interval in which stream is active
   Err += StreamConfig.get("UseStartEnd", NewStream->UseStartEnd);
   if (Err.isFail()) { // Start end flag not in config, assume false
      NewStream->UseStartEnd = false;
      Err.reset();
   }

   // Set Alarms for start and end time
   if (NewStream->UseStartEnd) {
      std::string StartTimeStr;
      std::string EndTimeStr;
      Err += StreamConfig.get("StartTime", StartTimeStr);
      CHECK_ERROR_ABORT(
          Err, "Stream {} requests UseStartEnd but no start time provided",
          StreamName);
      Err += StreamConfig.get("EndTime", EndTimeStr);
      CHECK_ERROR_ABORT(
          Err, "Stream {} requests UseStartEnd but no end time provided",
          StreamName);

      TimeInstant Start(StartTimeStr);
      TimeInstant End(EndTimeStr);
      std::string StartName = StreamName + "Start";
      std::string EndName   = StreamName + "End";
      NewStream->StartAlarm = Alarm(StartName, Start);
      NewStream->EndAlarm   = Alarm(EndName, End);
      ModelClock->attachAlarm(&(NewStream->StartAlarm));
      ModelClock->attachAlarm(&(NewStream->EndAlarm));
   } // endif UseStartEnd

   // Now we add the list of field names to the stream
   // First get the contents list
   std::vector<std::string> FieldContents;
   Err = StreamConfig.get("Contents", FieldContents);
   CHECK_ERROR_ABORT(Err, "Can not find contents for stream {}", StreamName);

   // If this is a write stream, add the time field for CF compliant
   // time information
   if (NewStream->Mode == IO::ModeWrite)
      NewStream->addField("time");

   // The contents are stored as an ordered set so we use the addField
   // interface to add each name. Note that in this context, the field
   // name can also be a group name. Group names are expanded during the
   // validate stage.
   for (int IField = 0; IField < FieldContents.size(); ++IField) {
      NewStream->addField(FieldContents[IField]);
   }

   // The contents list has not yet been validated.
   NewStream->Validated = false;

   // If we have made it to this point, we have a valid stream to add to
   // the list
   AllStreams[StreamName] = NewStream;

   return;

} // End IOStream create

//------------------------------------------------------------------------------
// Define all dimensions used. Returns a map of dimension names to defined
// dimension IDs.
void IOStream::defineAllDims(
    int FileID,                           ///< [in] id assigned to the IO file
    std::map<std::string, int> &AllDimIDs ///< [out] dim name, assigned ID
) {

   Error Err;

   for (auto IDim = Dimension::begin(); IDim != Dimension::end(); ++IDim) {
      std::string DimName = IDim->first;
      // For back compatibility, we also allow an older name (typically
      // with first name lower case). MaxCellsOnEdge is an exception since
      // it was named TWO in MPAS
      std::string OldDimName = DimName;
      if (DimName == "MaxCellsOnEdge") {
         OldDimName = "TWO";
      } else {
         OldDimName[0] = std::tolower(OldDimName[0]);
      }
      I4 Length = IDim->second->getLengthGlobal();
      I4 DimID;

      // First check to see if the dimension already exists in the file
      I4 InLength;
      Err = IO::getDimFromFile(FileID, DimName, DimID, InLength);
      if (Err.isFail()) { // dim not found
         // Try again using old name for back compatibility to MPAS
         Err = IO::getDimFromFile(FileID, OldDimName, DimID, InLength);
      }

      // If dim is found, use this DimID and check the length for consistency
      if (Err.isSuccess()) {
         if (InLength != Length and Length != IO::Unlimited)
            ABORT_ERROR("Inconsistent length for dimension {} in stream {}",
                        DimName, Name);

         // If dim is not found, define the dimension from the dim class if
         // it is a write operation, otherwise assume the dimension is not
         // needed
      } else {

         if (Mode == IO::Mode::ModeWrite) {
            DimID = IO::defineDim(FileID, DimName, Length);
         }

      } // end if found in file

      // Add the DimID to map for later use
      AllDimIDs[DimName] = DimID;

   } // end loop over all dims

   return;

} // End defineAllDims

//------------------------------------------------------------------------------
// Retrieves field size and dim lengths for non-distributed fields
// (distributed fields get this information from computeDecomp)
void IOStream::getFieldSize(
    std::shared_ptr<Field> FieldPtr, // [in] pointer to Field
    int &LocalSize,                  // [out] size of local array
    std::vector<int> &DimLengths     // [out] vector of local dim lengths
) {

   // Retrieve some basic field information
   std::string FieldName = FieldPtr->getName();
   int NDims             = FieldPtr->getNumDims();
   if (NDims == 0) { // scalar field
      LocalSize     = 1;
      DimLengths[0] = 1;
      return;
   } else if (NDims < 1) {
      ABORT_ERROR(
          "IOStream getFieldSize:Invalid number of dimensions for Field {}",
          FieldName);
   }

   std::vector<std::string> DimNames(NDims);
   FieldPtr->getDimNames(DimNames);

   LocalSize = 1;
   for (int IDim = 0; IDim < NDims; ++IDim) {
      std::string DimName                = DimNames[IDim];
      std::shared_ptr<Dimension> ThisDim = Dimension::get(DimName);
      DimLengths[IDim]                   = ThisDim->getLengthLocal();
      LocalSize *= DimLengths[IDim];
   }
}
//------------------------------------------------------------------------------
// Computes the parallel decomposition (offsets) for a field needed for parallel
// I/O. Return Decomp ID and array size for field.
void IOStream::computeDecomp(
    std::shared_ptr<Field> FieldPtr, // [in] pointer to Field
    int &DecompID,                   // [out] ID assigned to the decomposition
    int &LocalSize,                  // [out] size of local array
    std::vector<int> &DimLengths     // [out] vector of local dim lengths
) {

   // Retrieve some basic field information
   std::string FieldName   = FieldPtr->getName();
   IO::IODataType MyIOType = getFieldIOType(FieldPtr);
   int NDims               = FieldPtr->getNumDims();
   if (NDims < 0)
      ABORT_ERROR("IOStream computeDecomp: "
                  "Invalid number of dimensions for Field {}",
                  FieldName);

   std::vector<std::string> DimNames(NDims);
   FieldPtr->getDimNames(DimNames);

   // Get dimension and size information for each dimension
   // For the offset calculation, we also create dimension and dimension offsets
   // in index order up to the max of 5 dimensions. The active indices are
   // at the end in array index order.
   I4 GlobalSize        = 1;
   LocalSize            = 1;
   constexpr I4 MaxDims = 5;
   std::vector<I4> DimLengthsGlob(NDims, 1);
   std::vector<HostArray1DI4> DimOffsets(MaxDims);
   std::vector<I4> DimLengthGlob(MaxDims, 1); // lengths padded to MaxDims
   std::vector<I4> DimLengthLoc(MaxDims, 1);  // lengths padded to MaxDims

   for (int IDim = 0; IDim < NDims; ++IDim) {
      I4 StartDim                        = MaxDims - NDims;
      std::string DimName                = DimNames[IDim];
      std::shared_ptr<Dimension> ThisDim = Dimension::get(DimName);
      DimLengths[IDim]                   = ThisDim->getLengthLocal();
      DimLengthsGlob[IDim]               = ThisDim->getLengthGlobal();
      DimOffsets[StartDim + IDim]        = ThisDim->getOffset();
      DimLengthLoc[StartDim + IDim]      = DimLengths[IDim];
      DimLengthGlob[StartDim + IDim]     = DimLengthsGlob[IDim];
      LocalSize *= DimLengths[IDim];
      GlobalSize *= DimLengthsGlob[IDim];
   }

   // Create the data decomposition based on dimension information
   // Compute offset index (0-based global index of each element) in
   // linear address space. Needed for the decomposition definition.
   // -1 is used to denote indices that are not written
   // The linear offset along each dimension has already been computed
   // by the dimension class.

   std::vector<I4> Offset(LocalSize, -1);

   // Compute strides in linear space for each dimension
   // For the proper stride in storage order, we start from the rightmost
   // index in index order.
   std::vector<I4> StrideGlob(MaxDims, 1);
   std::vector<I4> StrideLoc(MaxDims, 1);
   for (int IDim = MaxDims - 2; IDim >= 0; --IDim) {
      StrideGlob[IDim] = StrideGlob[IDim + 1] * DimLengthGlob[IDim + 1];
      StrideLoc[IDim]  = StrideLoc[IDim + 1] * DimLengthLoc[IDim + 1];
   }

   // Compute full array offsets based on each dimensions linear offset
   // Skip -1 entries that are outside owned domain
   for (int N = 0; N < DimLengthLoc[0]; ++N) {
      I4 NGlob = 0;
      if (NDims >= 5)
         NGlob = (DimOffsets[0])(N);
      if (NGlob < 0)
         continue; // index outside owned domain
      for (int M = 0; M < DimLengthLoc[1]; ++M) {
         I4 MGlob = 0;
         if (NDims >= 4)
            MGlob = (DimOffsets[1])(M);
         if (MGlob < 0)
            continue;
         for (int K = 0; K < DimLengthLoc[2]; ++K) {
            I4 KGlob = 0;
            if (NDims >= 3)
               KGlob = (DimOffsets[2])(K);
            if (KGlob < 0)
               continue;
            for (int J = 0; J < DimLengthLoc[3]; ++J) {
               I4 JGlob = 0;
               if (NDims >= 2)
                  JGlob = (DimOffsets[3])(J);
               if (JGlob < 0)
                  continue;
               for (int I = 0; I < DimLengthLoc[4]; ++I) {
                  I4 IGlob = (DimOffsets[4])(I);
                  if (IGlob < 0)
                     continue;
                  int Add = I * StrideLoc[4] + J * StrideLoc[3] +
                            K * StrideLoc[2] + M * StrideLoc[1] +
                            N * StrideLoc[0];
                  Offset[Add] = IGlob * StrideGlob[4] + JGlob * StrideGlob[3] +
                                KGlob * StrideGlob[2] + MGlob * StrideGlob[1] +
                                NGlob * StrideGlob[0];
               }
            }
         }
      }
   }

   DecompID = IO::createDecomp(MyIOType, NDims, DimLengthsGlob, LocalSize,
                               Offset, IO::DefaultRearr);

   return;

} // End computeDecomp

//------------------------------------------------------------------------------
// Write all metadata associated with a field
void IOStream::writeFieldMeta(
    const std::string FieldName, // [in] metadata for field
    int FileID,                  // [in] id assigned to open file
    int FieldID                  // [in] id assigned to the field
) {

   // Get the field metadata entries
   std::shared_ptr<Field> ThisField  = Field::get(FieldName);
   std::shared_ptr<Metadata> AllMeta = ThisField->getAllMetadata();

   // Loop through all metadata - the Metadata type is an alias for a std::map
   // so we can use map iterators for the loop
   for (auto IMeta = AllMeta->begin(); IMeta != AllMeta->end(); ++IMeta) {

      // Get name
      std::string MetaName = IMeta->first;
      // Get value after determining the data type
      std::any MetaVal = IMeta->second;
      if (MetaVal.type() == typeid(I8)) {
         I8 MetaValI8 = std::any_cast<I8>(MetaVal);
         IO::writeMeta(MetaName, MetaValI8, FileID, FieldID);

      } else if (MetaVal.type() == typeid(I4)) {
         I4 MetaValI4 = std::any_cast<I4>(MetaVal);
         IO::writeMeta(MetaName, MetaValI4, FileID, FieldID);

      } else if (MetaVal.type() == typeid(R8)) {
         R8 MetaValR8 = std::any_cast<R8>(MetaVal);
         // if reduced precision is desired, convert to single (R4)
         if (ReducePrecision and !ThisField->retainPrecision()) {
            R4 MetaValR4 = MetaValR8;
            IO::writeMeta(MetaName, MetaValR4, FileID, FieldID);
         } else {
            IO::writeMeta(MetaName, MetaValR8, FileID, FieldID);
         }

      } else if (MetaVal.type() == typeid(R4)) {
         R4 MetaValR4 = std::any_cast<R4>(MetaVal);
         IO::writeMeta(MetaName, MetaValR4, FileID, FieldID);

      } else if (MetaVal.type() == typeid(bool)) {
         bool MetaValBool = std::any_cast<bool>(MetaVal);
         // bool is coerced to int in write
         IO::writeMeta(MetaName, MetaValBool, FileID, FieldID);

      } else if (MetaVal.type() == typeid(std::string)) {
         std::string MetaValStr = std::any_cast<std::string>(MetaVal);
         IO::writeMeta(MetaName, MetaValStr, FileID, FieldID);

         // If the metadata was assigned using a string literal, std::any
         // stores it as a char pointer so we need to convert it differently
      } else if (MetaVal.type() == typeid(const char *)) {
         const char *MetaValChar = std::any_cast<const char *>(MetaVal);
         IO::writeMeta(MetaName, MetaValChar, FileID, FieldID);

      } else { // unknown data type
         ABORT_ERROR("Unknown data type for Metadata {} in Field {}", MetaName,
                     FieldName);
      }

   } // end loop over metadata

   return;

} // End writeFieldMeta

//------------------------------------------------------------------------------
// Write a field's data array, performing any manipulations to reduce
// precision or move data between host and device
void IOStream::writeFieldData(
    std::shared_ptr<Field> FieldPtr,      // [in] field to write
    int FileID,                           // [in] id assigned to open file
    int FieldID,                          // [in] id assigned to the field
    std::map<std::string, int> &AllDimIDs // [in] dimension IDs
) {

   Error Err; // internal error code, default success

   // Retrieve some basic field information
   std::string FieldName = FieldPtr->getName();
   bool OnHost           = FieldPtr->isOnHost();
   bool IsDistributed    = FieldPtr->isDistributed();
   bool IsTimeDependent  = FieldPtr->isTimeDependent();
   bool RetainPrecision  = FieldPtr->retainPrecision();
   ArrayDataType MyType  = FieldPtr->getType();
   int NDims             = FieldPtr->getNumDims();
   if (NDims < 0)
      ABORT_ERROR("IOStream::writeFieldData: "
                  "Invalid number of dimensions for Field {}",
                  FieldName);

   // Create the decomposition needed for parallel I/O or if not decomposed
   // get the relevant size information
   int MyDecompID;
   int LocSize;
   int NDimsTmp = std::max(NDims, 1);
   std::vector<int> DimLengths(NDimsTmp);
   if (IsDistributed) {
      computeDecomp(FieldPtr, MyDecompID, LocSize, DimLengths);
   } else { // Get dimension lengths
      IOStream::getFieldSize(FieldPtr, LocSize, DimLengths);
      // Scalar data stored as an array with size 1 so reset the local NDims
      // to pick this up
      if (NDims == 0)
         ++NDims;
   }

   // Extract and write the array of data based on the type, dimension and
   // memory location. The IO routines require a contiguous data pointer on
   // the host. Kokkos array types do not guarantee contigous memory for
   // multi-dimensional arrays. Here we create a contiguous space and perform
   // any other transformations (host-device data transfer, reduce precision).
   void *DataPtr;
   void *FillValPtr;

   // Vector for contiguous storage - the appropriate vector will be
   // selected and resized later.
   std::vector<I4> DataI4(1);
   std::vector<I8> DataI8(1);
   std::vector<R4> DataR4(1);
   std::vector<R8> DataR8(1);
   I4 FillValI4;
   I8 FillValI8;
   R4 FillValR4;
   R8 FillValR8;

   switch (MyType) {

   // I4 Fields
   case ArrayDataType::I4:

      DataI4.resize(LocSize);
      DataPtr = DataI4.data();
      // get fill value
      Err += FieldPtr->getMetadata("FillValue", FillValI4);
      CHECK_ERROR_ABORT(Err, "Error retrieving FillValue for Field {}",
                        FieldName);
      FillValPtr = &FillValI4;

      switch (NDims) {
      case 1:
         if (OnHost) {
            HostArray1DI4 Data = FieldPtr->getDataArray<HostArray1DI4>();
            for (int I = 0; I < DimLengths[0]; ++I) {
               DataI4[I] = Data(I);
            }
         } else {
            Array1DI4 DataTmp  = FieldPtr->getDataArray<Array1DI4>();
            HostArray1DI4 Data = createHostMirrorCopy(DataTmp);
            for (int I = 0; I < DimLengths[0]; ++I) {
               DataI4[I] = Data(I);
            }
         }
         break;
      case 2:
         if (OnHost) {
            HostArray2DI4 Data = FieldPtr->getDataArray<HostArray2DI4>();
            int VecAdd         = 0;
            for (int J = 0; J < DimLengths[0]; ++J) {
               for (int I = 0; I < DimLengths[1]; ++I) {
                  DataI4[VecAdd] = Data(J, I);
                  ++VecAdd;
               }
            }
         } else {
            Array2DI4 DataTmp  = FieldPtr->getDataArray<Array2DI4>();
            HostArray2DI4 Data = createHostMirrorCopy(DataTmp);
            int VecAdd         = 0;
            for (int J = 0; J < DimLengths[0]; ++J) {
               for (int I = 0; I < DimLengths[1]; ++I) {
                  DataI4[VecAdd] = Data(J, I);
                  ++VecAdd;
               }
            }
         }
         break;
      case 3:
         if (OnHost) {
            HostArray3DI4 Data = FieldPtr->getDataArray<HostArray3DI4>();
            int VecAdd         = 0;
            for (int K = 0; K < DimLengths[0]; ++K) {
               for (int J = 0; J < DimLengths[1]; ++J) {
                  for (int I = 0; I < DimLengths[2]; ++I) {
                     DataI4[VecAdd] = Data(K, J, I);
                     ++VecAdd;
                  }
               }
            }
         } else {
            Array3DI4 DataTmp  = FieldPtr->getDataArray<Array3DI4>();
            HostArray3DI4 Data = createHostMirrorCopy(DataTmp);
            int VecAdd         = 0;
            for (int K = 0; K < DimLengths[0]; ++K) {
               for (int J = 0; J < DimLengths[1]; ++J) {
                  for (int I = 0; I < DimLengths[2]; ++I) {
                     DataI4[VecAdd] = Data(K, J, I);
                     ++VecAdd;
                  }
               }
            }
         }
         break;
      case 4:
         if (OnHost) {
            HostArray4DI4 Data = FieldPtr->getDataArray<HostArray4DI4>();
            int VecAdd         = 0;
            for (int L = 0; L < DimLengths[0]; ++L) {
               for (int K = 0; K < DimLengths[1]; ++K) {
                  for (int J = 0; J < DimLengths[2]; ++J) {
                     for (int I = 0; I < DimLengths[3]; ++I) {
                        DataI4[VecAdd] = Data(L, K, J, I);
                        ++VecAdd;
                     }
                  }
               }
            }
         } else {
            Array4DI4 DataTmp  = FieldPtr->getDataArray<Array4DI4>();
            HostArray4DI4 Data = createHostMirrorCopy(DataTmp);
            int VecAdd         = 0;
            for (int L = 0; L < DimLengths[0]; ++L) {
               for (int K = 0; K < DimLengths[1]; ++K) {
                  for (int J = 0; J < DimLengths[2]; ++J) {
                     for (int I = 0; I < DimLengths[3]; ++I) {
                        DataI4[VecAdd] = Data(L, K, J, I);
                        ++VecAdd;
                     }
                  }
               }
            }
         }
         break;
      case 5:
         if (OnHost) {
            HostArray5DI4 Data = FieldPtr->getDataArray<HostArray5DI4>();
            int VecAdd         = 0;
            for (int M = 0; M < DimLengths[0]; ++M) {
               for (int L = 0; L < DimLengths[1]; ++L) {
                  for (int K = 0; K < DimLengths[2]; ++K) {
                     for (int J = 0; J < DimLengths[3]; ++J) {
                        for (int I = 0; I < DimLengths[4]; ++I) {
                           DataI4[VecAdd] = Data(M, L, K, J, I);
                           ++VecAdd;
                        }
                     }
                  }
               }
            }
         } else {
            Array5DI4 DataTmp  = FieldPtr->getDataArray<Array5DI4>();
            HostArray5DI4 Data = createHostMirrorCopy(DataTmp);
            int VecAdd         = 0;
            for (int M = 0; M < DimLengths[0]; ++M) {
               for (int L = 0; L < DimLengths[1]; ++L) {
                  for (int K = 0; K < DimLengths[2]; ++K) {
                     for (int J = 0; J < DimLengths[3]; ++J) {
                        for (int I = 0; I < DimLengths[4]; ++I) {
                           DataI4[VecAdd] = Data(M, L, K, J, I);
                           ++VecAdd;
                        }
                     }
                  }
               }
            }
         }
         break;

      } // end switch NDims
      break; // end I4 type

   // I8 Fields
   case ArrayDataType::I8:

      DataI8.resize(LocSize);
      DataPtr = DataI8.data();
      // Get fill value
      Err += FieldPtr->getMetadata("FillValue", FillValI8);
      CHECK_ERROR_ABORT(Err, "Error retrieving FillValue for Field {}",
                        FieldName);
      FillValPtr = &FillValI8;

      switch (NDims) {
      case 1:
         if (OnHost) {
            HostArray1DI8 Data = FieldPtr->getDataArray<HostArray1DI8>();
            for (int I = 0; I < DimLengths[0]; ++I) {
               DataI8[I] = Data(I);
            }
         } else {
            Array1DI8 DataTmp  = FieldPtr->getDataArray<Array1DI8>();
            HostArray1DI8 Data = createHostMirrorCopy(DataTmp);
            for (int I = 0; I < DimLengths[0]; ++I) {
               DataI8[I] = Data(I);
            }
         }
         break;
      case 2:
         if (OnHost) {
            HostArray2DI8 Data = FieldPtr->getDataArray<HostArray2DI8>();
            int VecAdd         = 0;
            for (int J = 0; J < DimLengths[0]; ++J) {
               for (int I = 0; I < DimLengths[1]; ++I) {
                  DataI8[VecAdd] = Data(J, I);
                  ++VecAdd;
               }
            }
         } else {
            Array2DI8 DataTmp  = FieldPtr->getDataArray<Array2DI8>();
            HostArray2DI8 Data = createHostMirrorCopy(DataTmp);
            int VecAdd         = 0;
            for (int J = 0; J < DimLengths[0]; ++J) {
               for (int I = 0; I < DimLengths[1]; ++I) {
                  DataI8[VecAdd] = Data(J, I);
                  ++VecAdd;
               }
            }
         }
         break;
      case 3:
         if (OnHost) {
            HostArray3DI8 Data = FieldPtr->getDataArray<HostArray3DI8>();
            int VecAdd         = 0;
            for (int K = 0; K < DimLengths[0]; ++K) {
               for (int J = 0; J < DimLengths[1]; ++J) {
                  for (int I = 0; I < DimLengths[2]; ++I) {
                     DataI8[VecAdd] = Data(K, J, I);
                     ++VecAdd;
                  }
               }
            }
         } else {
            Array3DI8 DataTmp  = FieldPtr->getDataArray<Array3DI8>();
            HostArray3DI8 Data = createHostMirrorCopy(DataTmp);
            int VecAdd         = 0;
            for (int K = 0; K < DimLengths[0]; ++K) {
               for (int J = 0; J < DimLengths[1]; ++J) {
                  for (int I = 0; I < DimLengths[2]; ++I) {
                     DataI8[VecAdd] = Data(K, J, I);
                     ++VecAdd;
                  }
               }
            }
         }
         break;
      case 4:
         if (OnHost) {
            HostArray4DI8 Data = FieldPtr->getDataArray<HostArray4DI8>();
            int VecAdd         = 0;
            for (int L = 0; L < DimLengths[0]; ++L) {
               for (int K = 0; K < DimLengths[1]; ++K) {
                  for (int J = 0; J < DimLengths[2]; ++J) {
                     for (int I = 0; I < DimLengths[3]; ++I) {
                        DataI8[VecAdd] = Data(L, K, J, I);
                        ++VecAdd;
                     }
                  }
               }
            }
         } else {
            Array4DI8 DataTmp  = FieldPtr->getDataArray<Array4DI8>();
            HostArray4DI8 Data = createHostMirrorCopy(DataTmp);
            int VecAdd         = 0;
            for (int L = 0; L < DimLengths[0]; ++L) {
               for (int K = 0; K < DimLengths[1]; ++K) {
                  for (int J = 0; J < DimLengths[2]; ++J) {
                     for (int I = 0; I < DimLengths[3]; ++I) {
                        DataI8[VecAdd] = Data(L, K, J, I);
                        ++VecAdd;
                     }
                  }
               }
            }
         }
         break;
      case 5:
         if (OnHost) {
            HostArray5DI8 Data = FieldPtr->getDataArray<HostArray5DI8>();
            int VecAdd         = 0;
            for (int M = 0; M < DimLengths[0]; ++M) {
               for (int L = 0; L < DimLengths[1]; ++L) {
                  for (int K = 0; K < DimLengths[2]; ++K) {
                     for (int J = 0; J < DimLengths[3]; ++J) {
                        for (int I = 0; I < DimLengths[4]; ++I) {
                           DataI8[VecAdd] = Data(M, L, K, J, I);
                           ++VecAdd;
                        }
                     }
                  }
               }
            }
         } else {
            Array5DI8 DataTmp  = FieldPtr->getDataArray<Array5DI8>();
            HostArray5DI8 Data = createHostMirrorCopy(DataTmp);
            int VecAdd         = 0;
            for (int M = 0; M < DimLengths[0]; ++M) {
               for (int L = 0; L < DimLengths[1]; ++L) {
                  for (int K = 0; K < DimLengths[2]; ++K) {
                     for (int J = 0; J < DimLengths[3]; ++J) {
                        for (int I = 0; I < DimLengths[4]; ++I) {
                           DataI8[VecAdd] = Data(M, L, K, J, I);
                           ++VecAdd;
                        }
                     }
                  }
               }
            }
         }
         break;
      } // end switch NDims
      break; // end I8 type

   // R4 Fields
   case ArrayDataType::R4:

      DataR4.resize(LocSize);
      DataPtr = DataR4.data();
      // Get fill value
      Err += FieldPtr->getMetadata("FillValue", FillValR4);
      CHECK_ERROR_ABORT(Err, "Error retrieving FillValue for Field {}",
                        FieldName);
      FillValPtr = &FillValR4;

      switch (NDims) {
      case 1:
         if (OnHost) {
            HostArray1DR4 Data = FieldPtr->getDataArray<HostArray1DR4>();
            for (int I = 0; I < DimLengths[0]; ++I) {
               DataR4[I] = Data(I);
            }
         } else {
            Array1DR4 DataTmp  = FieldPtr->getDataArray<Array1DR4>();
            HostArray1DR4 Data = createHostMirrorCopy(DataTmp);
            for (int I = 0; I < DimLengths[0]; ++I) {
               DataR4[I] = Data(I);
            }
         }
         break;
      case 2:
         if (OnHost) {
            HostArray2DR4 Data = FieldPtr->getDataArray<HostArray2DR4>();
            int VecAdd         = 0;
            for (int J = 0; J < DimLengths[0]; ++J) {
               for (int I = 0; I < DimLengths[1]; ++I) {
                  DataR4[VecAdd] = Data(J, I);
                  ++VecAdd;
               }
            }
         } else {
            Array2DR4 DataTmp  = FieldPtr->getDataArray<Array2DR4>();
            HostArray2DR4 Data = createHostMirrorCopy(DataTmp);
            int VecAdd         = 0;
            for (int J = 0; J < DimLengths[0]; ++J) {
               for (int I = 0; I < DimLengths[1]; ++I) {
                  DataR4[VecAdd] = Data(J, I);
                  ++VecAdd;
               }
            }
         }
         break;
      case 3:
         if (OnHost) {
            HostArray3DR4 Data = FieldPtr->getDataArray<HostArray3DR4>();
            int VecAdd         = 0;
            for (int K = 0; K < DimLengths[0]; ++K) {
               for (int J = 0; J < DimLengths[1]; ++J) {
                  for (int I = 0; I < DimLengths[2]; ++I) {
                     DataR4[VecAdd] = Data(K, J, I);
                     ++VecAdd;
                  }
               }
            }
         } else {
            Array3DR4 DataTmp  = FieldPtr->getDataArray<Array3DR4>();
            HostArray3DR4 Data = createHostMirrorCopy(DataTmp);
            int VecAdd         = 0;
            for (int K = 0; K < DimLengths[0]; ++K) {
               for (int J = 0; J < DimLengths[1]; ++J) {
                  for (int I = 0; I < DimLengths[2]; ++I) {
                     DataR4[VecAdd] = Data(K, J, I);
                     ++VecAdd;
                  }
               }
            }
         }
         break;
      case 4:
         if (OnHost) {
            HostArray4DR4 Data = FieldPtr->getDataArray<HostArray4DR4>();
            int VecAdd         = 0;
            for (int L = 0; L < DimLengths[0]; ++L) {
               for (int K = 0; K < DimLengths[1]; ++K) {
                  for (int J = 0; J < DimLengths[2]; ++J) {
                     for (int I = 0; I < DimLengths[3]; ++I) {
                        DataR4[VecAdd] = Data(L, K, J, I);
                        ++VecAdd;
                     }
                  }
               }
            }
         } else {
            Array4DR4 DataTmp  = FieldPtr->getDataArray<Array4DR4>();
            HostArray4DR4 Data = createHostMirrorCopy(DataTmp);
            int VecAdd         = 0;
            for (int L = 0; L < DimLengths[0]; ++L) {
               for (int K = 0; K < DimLengths[1]; ++K) {
                  for (int J = 0; J < DimLengths[2]; ++J) {
                     for (int I = 0; I < DimLengths[3]; ++I) {
                        DataR4[VecAdd] = Data(L, K, J, I);
                        ++VecAdd;
                     }
                  }
               }
            }
         }
         break;
      case 5:
         if (OnHost) {
            HostArray5DR4 Data = FieldPtr->getDataArray<HostArray5DR4>();
            int VecAdd         = 0;
            for (int M = 0; M < DimLengths[0]; ++M) {
               for (int L = 0; L < DimLengths[1]; ++L) {
                  for (int K = 0; K < DimLengths[2]; ++K) {
                     for (int J = 0; J < DimLengths[3]; ++J) {
                        for (int I = 0; I < DimLengths[4]; ++I) {
                           DataR4[VecAdd] = Data(M, L, K, J, I);
                           ++VecAdd;
                        }
                     }
                  }
               }
            }
         } else {
            Array5DR4 DataTmp  = FieldPtr->getDataArray<Array5DR4>();
            HostArray5DR4 Data = createHostMirrorCopy(DataTmp);
            int VecAdd         = 0;
            for (int M = 0; M < DimLengths[0]; ++M) {
               for (int L = 0; L < DimLengths[1]; ++L) {
                  for (int K = 0; K < DimLengths[2]; ++K) {
                     for (int J = 0; J < DimLengths[3]; ++J) {
                        for (int I = 0; I < DimLengths[4]; ++I) {
                           DataR4[VecAdd] = Data(M, L, K, J, I);
                           ++VecAdd;
                        }
                     }
                  }
               }
            }
         }
         break;
      } // end switch NDims
      break; // end R4 type

   // R8 Fields
   case ArrayDataType::R8:

      // Get fill value
      Err += FieldPtr->getMetadata("FillValue", FillValR8);
      CHECK_ERROR_ABORT(Err, "Error retrieving FillValue for Field {}",
                        FieldName);
      DataR8.resize(LocSize);
      if (ReducePrecision and !RetainPrecision) {
         FillValR4  = FillValR8;
         FillValPtr = &FillValR4;
         DataR4.resize(LocSize);
         DataPtr = DataR4.data();
      } else {
         FillValPtr = &FillValR8;
         DataPtr    = DataR8.data();
      }

      switch (NDims) {
      case 1:
         if (OnHost) {
            HostArray1DR8 Data = FieldPtr->getDataArray<HostArray1DR8>();
            for (int I = 0; I < DimLengths[0]; ++I) {
               DataR8[I] = Data(I);
            }
         } else {
            Array1DR8 DataTmp  = FieldPtr->getDataArray<Array1DR8>();
            HostArray1DR8 Data = createHostMirrorCopy(DataTmp);
            for (int I = 0; I < DimLengths[0]; ++I) {
               DataR8[I] = Data(I);
            }
         }
         break;
      case 2:
         if (OnHost) {
            HostArray2DR8 Data = FieldPtr->getDataArray<HostArray2DR8>();
            int VecAdd         = 0;
            for (int J = 0; J < DimLengths[0]; ++J) {
               for (int I = 0; I < DimLengths[1]; ++I) {
                  DataR8[VecAdd] = Data(J, I);
                  ++VecAdd;
               }
            }
         } else {
            Array2DR8 DataTmp  = FieldPtr->getDataArray<Array2DR8>();
            HostArray2DR8 Data = createHostMirrorCopy(DataTmp);
            int VecAdd         = 0;
            for (int J = 0; J < DimLengths[0]; ++J) {
               for (int I = 0; I < DimLengths[1]; ++I) {
                  DataR8[VecAdd] = Data(J, I);
                  ++VecAdd;
               }
            }
         }
         break;
      case 3:
         if (OnHost) {
            HostArray3DR8 Data = FieldPtr->getDataArray<HostArray3DR8>();
            int VecAdd         = 0;
            for (int K = 0; K < DimLengths[0]; ++K) {
               for (int J = 0; J < DimLengths[1]; ++J) {
                  for (int I = 0; I < DimLengths[2]; ++I) {
                     DataR8[VecAdd] = Data(K, J, I);
                     ++VecAdd;
                  }
               }
            }
         } else {
            Array3DR8 DataTmp  = FieldPtr->getDataArray<Array3DR8>();
            HostArray3DR8 Data = createHostMirrorCopy(DataTmp);
            int VecAdd         = 0;
            for (int K = 0; K < DimLengths[0]; ++K) {
               for (int J = 0; J < DimLengths[1]; ++J) {
                  for (int I = 0; I < DimLengths[2]; ++I) {
                     DataR8[VecAdd] = Data(K, J, I);
                     ++VecAdd;
                  }
               }
            }
         }
         break;
      case 4:
         if (OnHost) {
            HostArray4DR8 Data = FieldPtr->getDataArray<HostArray4DR8>();
            int VecAdd         = 0;
            for (int L = 0; L < DimLengths[0]; ++L) {
               for (int K = 0; K < DimLengths[1]; ++K) {
                  for (int J = 0; J < DimLengths[2]; ++J) {
                     for (int I = 0; I < DimLengths[3]; ++I) {
                        DataR8[VecAdd] = Data(L, K, J, I);
                        ++VecAdd;
                     }
                  }
               }
            }
         } else {
            Array4DR8 DataTmp  = FieldPtr->getDataArray<Array4DR8>();
            HostArray4DR8 Data = createHostMirrorCopy(DataTmp);
            int VecAdd         = 0;
            for (int L = 0; L < DimLengths[0]; ++L) {
               for (int K = 0; K < DimLengths[1]; ++K) {
                  for (int J = 0; J < DimLengths[2]; ++J) {
                     for (int I = 0; I < DimLengths[3]; ++I) {
                        DataR8[VecAdd] = Data(L, K, J, I);
                        ++VecAdd;
                     }
                  }
               }
            }
         }
         break;
      case 5:
         if (OnHost) {
            HostArray5DR8 Data = FieldPtr->getDataArray<HostArray5DR8>();
            int VecAdd         = 0;
            for (int M = 0; M < DimLengths[0]; ++M) {
               for (int L = 0; L < DimLengths[1]; ++L) {
                  for (int K = 0; K < DimLengths[2]; ++K) {
                     for (int J = 0; J < DimLengths[3]; ++J) {
                        for (int I = 0; I < DimLengths[4]; ++I) {
                           DataR8[VecAdd] = Data(M, L, K, J, I);
                           ++VecAdd;
                        }
                     }
                  }
               }
            }
         } else {
            Array5DR8 DataTmp  = FieldPtr->getDataArray<Array5DR8>();
            HostArray5DR8 Data = createHostMirrorCopy(DataTmp);
            int VecAdd         = 0;
            for (int M = 0; M < DimLengths[0]; ++M) {
               for (int L = 0; L < DimLengths[1]; ++L) {
                  for (int K = 0; K < DimLengths[2]; ++K) {
                     for (int J = 0; J < DimLengths[3]; ++J) {
                        for (int I = 0; I < DimLengths[4]; ++I) {
                           DataR8[VecAdd] = Data(M, L, K, J, I);
                           ++VecAdd;
                        }
                     }
                  }
               }
            }
         }
         break;
      } // end switch NDims
      if (ReducePrecision and !RetainPrecision) {
         for (int I = 0; I < LocSize; ++I) {
            DataR4[I] = DataR8[I];
         }
      }
      break; // end R8 type

   default:
      ABORT_ERROR("Cannot determine data type for field {}", FieldName);

   } // end switch data type

   // In the case where the field is not time-dependent but this is
   // a multi-slice file, we reset the Frame temporarily for this field
   int FldFrame = Frame; // default time slice for file
   if (!IsTimeDependent)
      FldFrame = -1;

   // Write the data
   if (IsDistributed) {
      IO::writeArray(DataPtr, LocSize, FillValPtr, FileID, MyDecompID, FieldID,
                     FldFrame);

      // Clean up the decomp
      IO::destroyDecomp(MyDecompID);

   } else {
      IO::writeNDVar(DataPtr, FileID, FieldID, FldFrame, &DimLengths);
   }

   return;

} // end writeFieldData

//------------------------------------------------------------------------------
// Read a field's data array, performing any manipulations to reduce
// precision or move data between host and device
Error IOStream::readFieldData(
    std::shared_ptr<Field> FieldPtr,       // [in] field to read
    int FileID,                            // [in] id assigned to open file
    std::map<std::string, int> &AllDimIDs, // [in] dimension IDs
    int &FieldID                           // [out] id assigned to the field
) {

   Error Err; // Error code for IO routines

   // Retrieve some basic field information
   std::string FieldName = FieldPtr->getName();
   // For MPAS back compatibility, the old name has a first letter that is
   // lower case
   std::string OldFieldName = FieldName;
   OldFieldName[0]          = std::tolower(OldFieldName[0]);
   // Check for "Layer" in field name, backwards compatibility requires
   // replacing "Layer" with "Level"
   std::string OmegaSubStr = "Layer";
   std::string MPASSubStr  = "Level";
   size_t pos              = OldFieldName.find(OmegaSubStr);
   if (pos != std::string::npos) {
      OldFieldName.replace(pos, OmegaSubStr.length(), MPASSubStr);
   }
   bool OnHost          = FieldPtr->isOnHost();
   bool IsDistributed   = FieldPtr->isDistributed();
   bool IsTimeDependent = FieldPtr->isTimeDependent();
   ArrayDataType MyType = FieldPtr->getType();
   int NDims            = FieldPtr->getNumDims();
   if (NDims < 0)
      ABORT_ERROR("IOStream readFieldData: "
                  "Invalid number of dimensions for Field {}",
                  FieldName);

   // Create the decomposition needed for parallel I/O or if not decomposed
   // get the relevant size information
   int DecompID;
   int LocSize;
   int NDimsTmp = std::max(NDims, 1);
   std::vector<int> DimLengths(NDimsTmp);
   if (IsDistributed) {
      computeDecomp(FieldPtr, DecompID, LocSize, DimLengths);
   } else { // Get dimension lengths
      IOStream::getFieldSize(FieldPtr, LocSize, DimLengths);
      // Scalar data stored as an array with size 1 so reset the local NDims
      // to pick this up
      if (NDims == 0)
         ++NDims;
   }

   // The IO routines require a pointer to a contiguous memory on the host
   // so we first read into a vector. Only one of the vectors below will
   // be used and resized appropriately.
   void *DataPtr;
   std::vector<I4> DataI4(1);
   std::vector<I8> DataI8(1);
   std::vector<R4> DataR4(1);
   std::vector<R8> DataR8(1);

   switch (MyType) {
   case ArrayDataType::I4:
      DataI4.resize(LocSize);
      DataPtr = DataI4.data();
      break;
   case ArrayDataType::I8:
      DataI8.resize(LocSize);
      DataPtr = DataI8.data();
      break;
   case ArrayDataType::R4:
      DataR4.resize(LocSize);
      DataPtr = DataR4.data();
      break;
   case ArrayDataType::R8:
      DataR8.resize(LocSize);
      DataPtr = DataR8.data();
      break;
   case ArrayDataType::Unknown:
      ABORT_ERROR("IOStream readFieldData: Unknown data array type");
   default:
      ABORT_ERROR("IOStream readFieldData: Invalid data array type");
   }

   // read data into vector
   if (IsDistributed) {
      Err =
          IO::readArray(DataPtr, LocSize, FieldName, FileID, DecompID, FieldID);
   } else {
      Err = IO::readNDVar(DataPtr, FieldName, FileID, FieldID);
   }
   // For back compatibility, try to read again with old field name
   if (Err.isFail()) {
      if (IsDistributed) {
         Err = IO::readArray(DataPtr, LocSize, OldFieldName, FileID, DecompID,
                             FieldID);
      } else {
         Err = IO::readNDVar(DataPtr, OldFieldName, FileID, FieldID);
      }
      if (Err.isFail()) // still cannot find field, return with error
         RETURN_ERROR(
             Err, ErrorCode::Fail,
             "IOStream::readFieldData: Field {} or {} not found in stream {}",
             FieldName, OldFieldName, Name);
   }

   // Unpack vector into array based on type, dims and location
   switch (MyType) {

   // I4 Fields
   case ArrayDataType::I4:
      switch (NDims) {
      case 1:
         if (OnHost) {
            HostArray1DI4 Data = FieldPtr->getDataArray<HostArray1DI4>();
            for (int I = 0; I < DimLengths[0]; ++I) {
               Data(I) = DataI4[I];
            }
         } else {
            Array1DI4 DataTmp  = FieldPtr->getDataArray<Array1DI4>();
            HostArray1DI4 Data = createHostMirrorCopy(DataTmp);
            for (int I = 0; I < DimLengths[0]; ++I) {
               Data(I) = DataI4[I];
            }
            deepCopy(DataTmp, Data);
         }
         break;
      case 2:
         if (OnHost) {
            HostArray2DI4 Data = FieldPtr->getDataArray<HostArray2DI4>();
            int VecAdd         = 0;
            for (int J = 0; J < DimLengths[0]; ++J) {
               for (int I = 0; I < DimLengths[1]; ++I) {
                  Data(J, I) = DataI4[VecAdd];
                  ++VecAdd;
               }
            }
         } else {
            Array2DI4 DataTmp  = FieldPtr->getDataArray<Array2DI4>();
            HostArray2DI4 Data = createHostMirrorCopy(DataTmp);
            int VecAdd         = 0;
            for (int J = 0; J < DimLengths[0]; ++J) {
               for (int I = 0; I < DimLengths[1]; ++I) {
                  Data(J, I) = DataI4[VecAdd];
                  ++VecAdd;
               }
            }
            deepCopy(DataTmp, Data);
         }
         break;
      case 3:
         if (OnHost) {
            HostArray3DI4 Data = FieldPtr->getDataArray<HostArray3DI4>();
            int VecAdd         = 0;
            for (int K = 0; K < DimLengths[0]; ++K) {
               for (int J = 0; J < DimLengths[1]; ++J) {
                  for (int I = 0; I < DimLengths[2]; ++I) {
                     Data(K, J, I) = DataI4[VecAdd];
                     ++VecAdd;
                  }
               }
            }
         } else {
            Array3DI4 DataTmp  = FieldPtr->getDataArray<Array3DI4>();
            HostArray3DI4 Data = createHostMirrorCopy(DataTmp);
            int VecAdd         = 0;
            for (int K = 0; K < DimLengths[0]; ++K) {
               for (int J = 0; J < DimLengths[1]; ++J) {
                  for (int I = 0; I < DimLengths[2]; ++I) {
                     Data(K, J, I) = DataI4[VecAdd];
                     ++VecAdd;
                  }
               }
            }
            deepCopy(DataTmp, Data);
         }
         break;
      case 4:
         if (OnHost) {
            HostArray4DI4 Data = FieldPtr->getDataArray<HostArray4DI4>();
            int VecAdd         = 0;
            for (int L = 0; L < DimLengths[0]; ++L) {
               for (int K = 0; K < DimLengths[1]; ++K) {
                  for (int J = 0; J < DimLengths[2]; ++J) {
                     for (int I = 0; I < DimLengths[3]; ++I) {
                        Data(L, K, J, I) = DataI4[VecAdd];
                        ++VecAdd;
                     }
                  }
               }
            }
         } else {
            Array4DI4 DataTmp  = FieldPtr->getDataArray<Array4DI4>();
            HostArray4DI4 Data = createHostMirrorCopy(DataTmp);
            int VecAdd         = 0;
            for (int L = 0; L < DimLengths[0]; ++L) {
               for (int K = 0; K < DimLengths[1]; ++K) {
                  for (int J = 0; J < DimLengths[2]; ++J) {
                     for (int I = 0; I < DimLengths[3]; ++I) {
                        Data(L, K, J, I) = DataI4[VecAdd];
                        ++VecAdd;
                     }
                  }
               }
            }
            deepCopy(DataTmp, Data);
         }
         break;
      case 5:
         if (OnHost) {
            HostArray5DI4 Data = FieldPtr->getDataArray<HostArray5DI4>();
            int VecAdd         = 0;
            for (int M = 0; M < DimLengths[0]; ++M) {
               for (int L = 0; L < DimLengths[1]; ++L) {
                  for (int K = 0; K < DimLengths[2]; ++K) {
                     for (int J = 0; J < DimLengths[3]; ++J) {
                        for (int I = 0; I < DimLengths[4]; ++I) {
                           Data(M, L, K, J, I) = DataI4[VecAdd];
                           ++VecAdd;
                        }
                     }
                  }
               }
            }
         } else {
            Array5DI4 DataTmp  = FieldPtr->getDataArray<Array5DI4>();
            HostArray5DI4 Data = createHostMirrorCopy(DataTmp);
            int VecAdd         = 0;
            for (int M = 0; M < DimLengths[0]; ++M) {
               for (int L = 0; L < DimLengths[1]; ++L) {
                  for (int K = 0; K < DimLengths[2]; ++K) {
                     for (int J = 0; J < DimLengths[3]; ++J) {
                        for (int I = 0; I < DimLengths[4]; ++I) {
                           Data(M, L, K, J, I) = DataI4[VecAdd];
                           ++VecAdd;
                        }
                     }
                  }
               }
            }
            deepCopy(DataTmp, Data);
         }
         break;
      } // end switch NDims
      break; // end I4 fields

   // I8 Fields
   case ArrayDataType::I8:
      switch (NDims) {
      case 1:
         if (OnHost) {
            HostArray1DI8 Data = FieldPtr->getDataArray<HostArray1DI8>();
            for (int I = 0; I < DimLengths[0]; ++I) {
               Data(I) = DataI8[I];
            }
         } else {
            Array1DI8 DataTmp  = FieldPtr->getDataArray<Array1DI8>();
            HostArray1DI8 Data = createHostMirrorCopy(DataTmp);
            for (int I = 0; I < DimLengths[0]; ++I) {
               Data(I) = DataI8[I];
            }
            deepCopy(DataTmp, Data);
         }
         break;
      case 2:
         if (OnHost) {
            HostArray2DI8 Data = FieldPtr->getDataArray<HostArray2DI8>();
            int VecAdd         = 0;
            for (int J = 0; J < DimLengths[0]; ++J) {
               for (int I = 0; I < DimLengths[1]; ++I) {
                  Data(J, I) = DataI8[VecAdd];
                  ++VecAdd;
               }
            }
         } else {
            Array2DI8 DataTmp  = FieldPtr->getDataArray<Array2DI8>();
            HostArray2DI8 Data = createHostMirrorCopy(DataTmp);
            int VecAdd         = 0;
            for (int J = 0; J < DimLengths[0]; ++J) {
               for (int I = 0; I < DimLengths[1]; ++I) {
                  Data(J, I) = DataI8[VecAdd];
                  ++VecAdd;
               }
            }
            deepCopy(DataTmp, Data);
         }
         break;
      case 3:
         if (OnHost) {
            HostArray3DI8 Data = FieldPtr->getDataArray<HostArray3DI8>();
            int VecAdd         = 0;
            for (int K = 0; K < DimLengths[0]; ++K) {
               for (int J = 0; J < DimLengths[1]; ++J) {
                  for (int I = 0; I < DimLengths[2]; ++I) {
                     Data(K, J, I) = DataI8[VecAdd];
                     ++VecAdd;
                  }
               }
            }
         } else {
            Array3DI8 DataTmp  = FieldPtr->getDataArray<Array3DI8>();
            HostArray3DI8 Data = createHostMirrorCopy(DataTmp);
            int VecAdd         = 0;
            for (int K = 0; K < DimLengths[0]; ++K) {
               for (int J = 0; J < DimLengths[1]; ++J) {
                  for (int I = 0; I < DimLengths[2]; ++I) {
                     Data(K, J, I) = DataI8[VecAdd];
                     ++VecAdd;
                  }
               }
            }
            deepCopy(DataTmp, Data);
         }
         break;
      case 4:
         if (OnHost) {
            HostArray4DI8 Data = FieldPtr->getDataArray<HostArray4DI8>();
            int VecAdd         = 0;
            for (int L = 0; L < DimLengths[0]; ++L) {
               for (int K = 0; K < DimLengths[1]; ++K) {
                  for (int J = 0; J < DimLengths[2]; ++J) {
                     for (int I = 0; I < DimLengths[3]; ++I) {
                        Data(L, K, J, I) = DataI8[VecAdd];
                        ++VecAdd;
                     }
                  }
               }
            }
         } else {
            Array4DI8 DataTmp  = FieldPtr->getDataArray<Array4DI8>();
            HostArray4DI8 Data = createHostMirrorCopy(DataTmp);
            int VecAdd         = 0;
            for (int L = 0; L < DimLengths[0]; ++L) {
               for (int K = 0; K < DimLengths[1]; ++K) {
                  for (int J = 0; J < DimLengths[2]; ++J) {
                     for (int I = 0; I < DimLengths[3]; ++I) {
                        Data(L, K, J, I) = DataI8[VecAdd];
                        ++VecAdd;
                     }
                  }
               }
            }
            deepCopy(DataTmp, Data);
         }
         break;
      case 5:
         if (OnHost) {
            HostArray5DI8 Data = FieldPtr->getDataArray<HostArray5DI8>();
            int VecAdd         = 0;
            for (int M = 0; M < DimLengths[0]; ++M) {
               for (int L = 0; L < DimLengths[1]; ++L) {
                  for (int K = 0; K < DimLengths[2]; ++K) {
                     for (int J = 0; J < DimLengths[3]; ++J) {
                        for (int I = 0; I < DimLengths[4]; ++I) {
                           Data(M, L, K, J, I) = DataI8[VecAdd];
                           ++VecAdd;
                        }
                     }
                  }
               }
            }
         } else {
            Array5DI8 DataTmp  = FieldPtr->getDataArray<Array5DI8>();
            HostArray5DI8 Data = createHostMirrorCopy(DataTmp);
            int VecAdd         = 0;
            for (int M = 0; M < DimLengths[0]; ++M) {
               for (int L = 0; L < DimLengths[1]; ++L) {
                  for (int K = 0; K < DimLengths[2]; ++K) {
                     for (int J = 0; J < DimLengths[3]; ++J) {
                        for (int I = 0; I < DimLengths[4]; ++I) {
                           Data(M, L, K, J, I) = DataI8[VecAdd];
                           ++VecAdd;
                        }
                     }
                  }
               }
            }
            deepCopy(DataTmp, Data);
         }
         break;
      } // end switch NDims
      break; // end I8 fields

   // R4 Fields
   case ArrayDataType::R4:
      switch (NDims) {
      case 1:
         if (OnHost) {
            HostArray1DR4 Data = FieldPtr->getDataArray<HostArray1DR4>();
            for (int I = 0; I < DimLengths[0]; ++I) {
               Data(I) = DataR4[I];
            }
         } else {
            Array1DR4 DataTmp  = FieldPtr->getDataArray<Array1DR4>();
            HostArray1DR4 Data = createHostMirrorCopy(DataTmp);
            for (int I = 0; I < DimLengths[0]; ++I) {
               Data(I) = DataR4[I];
            }
            deepCopy(DataTmp, Data);
         }
         break;
      case 2:
         if (OnHost) {
            HostArray2DR4 Data = FieldPtr->getDataArray<HostArray2DR4>();
            int VecAdd         = 0;
            for (int J = 0; J < DimLengths[0]; ++J) {
               for (int I = 0; I < DimLengths[1]; ++I) {
                  Data(J, I) = DataR4[VecAdd];
                  ++VecAdd;
               }
            }
         } else {
            Array2DR4 DataTmp  = FieldPtr->getDataArray<Array2DR4>();
            HostArray2DR4 Data = createHostMirrorCopy(DataTmp);
            int VecAdd         = 0;
            for (int J = 0; J < DimLengths[0]; ++J) {
               for (int I = 0; I < DimLengths[1]; ++I) {
                  Data(J, I) = DataR4[VecAdd];
                  ++VecAdd;
               }
            }
            deepCopy(DataTmp, Data);
         }
         break;
      case 3:
         if (OnHost) {
            HostArray3DR4 Data = FieldPtr->getDataArray<HostArray3DR4>();
            int VecAdd         = 0;
            for (int K = 0; K < DimLengths[0]; ++K) {
               for (int J = 0; J < DimLengths[1]; ++J) {
                  for (int I = 0; I < DimLengths[2]; ++I) {
                     Data(K, J, I) = DataR4[VecAdd];
                     ++VecAdd;
                  }
               }
            }
         } else {
            Array3DR4 DataTmp  = FieldPtr->getDataArray<Array3DR4>();
            HostArray3DR4 Data = createHostMirrorCopy(DataTmp);
            int VecAdd         = 0;
            for (int K = 0; K < DimLengths[0]; ++K) {
               for (int J = 0; J < DimLengths[1]; ++J) {
                  for (int I = 0; I < DimLengths[2]; ++I) {
                     Data(K, J, I) = DataR4[VecAdd];
                     ++VecAdd;
                  }
               }
            }
            deepCopy(DataTmp, Data);
         }
         break;
      case 4:
         if (OnHost) {
            HostArray4DR4 Data = FieldPtr->getDataArray<HostArray4DR4>();
            int VecAdd         = 0;
            for (int L = 0; L < DimLengths[0]; ++L) {
               for (int K = 0; K < DimLengths[1]; ++K) {
                  for (int J = 0; J < DimLengths[2]; ++J) {
                     for (int I = 0; I < DimLengths[3]; ++I) {
                        Data(L, K, J, I) = DataR4[VecAdd];
                        ++VecAdd;
                     }
                  }
               }
            }
         } else {
            Array4DR4 DataTmp  = FieldPtr->getDataArray<Array4DR4>();
            HostArray4DR4 Data = createHostMirrorCopy(DataTmp);
            int VecAdd         = 0;
            for (int L = 0; L < DimLengths[0]; ++L) {
               for (int K = 0; K < DimLengths[1]; ++K) {
                  for (int J = 0; J < DimLengths[2]; ++J) {
                     for (int I = 0; I < DimLengths[3]; ++I) {
                        Data(L, K, J, I) = DataR4[VecAdd];
                        ++VecAdd;
                     }
                  }
               }
            }
            deepCopy(DataTmp, Data);
         }
         break;
      case 5:
         if (OnHost) {
            HostArray5DR4 Data = FieldPtr->getDataArray<HostArray5DR4>();
            int VecAdd         = 0;
            for (int M = 0; M < DimLengths[0]; ++M) {
               for (int L = 0; L < DimLengths[1]; ++L) {
                  for (int K = 0; K < DimLengths[2]; ++K) {
                     for (int J = 0; J < DimLengths[3]; ++J) {
                        for (int I = 0; I < DimLengths[4]; ++I) {
                           Data(M, L, K, J, I) = DataR4[VecAdd];
                           ++VecAdd;
                        }
                     }
                  }
               }
            }
         } else {
            Array5DR4 DataTmp  = FieldPtr->getDataArray<Array5DR4>();
            HostArray5DR4 Data = createHostMirrorCopy(DataTmp);
            int VecAdd         = 0;
            for (int M = 0; M < DimLengths[0]; ++M) {
               for (int L = 0; L < DimLengths[1]; ++L) {
                  for (int K = 0; K < DimLengths[2]; ++K) {
                     for (int J = 0; J < DimLengths[3]; ++J) {
                        for (int I = 0; I < DimLengths[4]; ++I) {
                           Data(M, L, K, J, I) = DataR4[VecAdd];
                           ++VecAdd;
                        }
                     }
                  }
               }
            }
            deepCopy(DataTmp, Data);
         }
         break;
      } // end switch NDims
      break; // end R4 fields

   // R8 Fields
   case ArrayDataType::R8:
      switch (NDims) {
      case 1:
         if (OnHost) {
            HostArray1DR8 Data = FieldPtr->getDataArray<HostArray1DR8>();
            for (int I = 0; I < DimLengths[0]; ++I) {
               Data(I) = DataR8[I];
            }
         } else {
            Array1DR8 DataTmp  = FieldPtr->getDataArray<Array1DR8>();
            HostArray1DR8 Data = createHostMirrorCopy(DataTmp);
            for (int I = 0; I < DimLengths[0]; ++I) {
               Data(I) = DataR8[I];
            }
            deepCopy(DataTmp, Data);
         }
         break;
      case 2:
         if (OnHost) {
            HostArray2DR8 Data = FieldPtr->getDataArray<HostArray2DR8>();
            int VecAdd         = 0;
            for (int J = 0; J < DimLengths[0]; ++J) {
               for (int I = 0; I < DimLengths[1]; ++I) {
                  Data(J, I) = DataR8[VecAdd];
                  ++VecAdd;
               }
            }
         } else {
            Array2DR8 DataTmp  = FieldPtr->getDataArray<Array2DR8>();
            HostArray2DR8 Data = createHostMirrorCopy(DataTmp);
            int VecAdd         = 0;
            for (int J = 0; J < DimLengths[0]; ++J) {
               for (int I = 0; I < DimLengths[1]; ++I) {
                  Data(J, I) = DataR8[VecAdd];
                  ++VecAdd;
               }
            }
            deepCopy(DataTmp, Data);
         }
         break;
      case 3:
         if (OnHost) {
            HostArray3DR8 Data = FieldPtr->getDataArray<HostArray3DR8>();
            int VecAdd         = 0;
            for (int K = 0; K < DimLengths[0]; ++K) {
               for (int J = 0; J < DimLengths[1]; ++J) {
                  for (int I = 0; I < DimLengths[2]; ++I) {
                     Data(K, J, I) = DataR8[VecAdd];
                     ++VecAdd;
                  }
               }
            }
         } else {
            Array3DR8 DataTmp  = FieldPtr->getDataArray<Array3DR8>();
            HostArray3DR8 Data = createHostMirrorCopy(DataTmp);
            int VecAdd         = 0;
            for (int K = 0; K < DimLengths[0]; ++K) {
               for (int J = 0; J < DimLengths[1]; ++J) {
                  for (int I = 0; I < DimLengths[2]; ++I) {
                     Data(K, J, I) = DataR8[VecAdd];
                     ++VecAdd;
                  }
               }
            }
            deepCopy(DataTmp, Data);
         }
         break;
      case 4:
         if (OnHost) {
            HostArray4DR8 Data = FieldPtr->getDataArray<HostArray4DR8>();
            int VecAdd         = 0;
            for (int L = 0; L < DimLengths[0]; ++L) {
               for (int K = 0; K < DimLengths[1]; ++K) {
                  for (int J = 0; J < DimLengths[2]; ++J) {
                     for (int I = 0; I < DimLengths[3]; ++I) {
                        Data(L, K, J, I) = DataR8[VecAdd];
                        ++VecAdd;
                     }
                  }
               }
            }
         } else {
            Array4DR8 DataTmp  = FieldPtr->getDataArray<Array4DR8>();
            HostArray4DR8 Data = createHostMirrorCopy(DataTmp);
            int VecAdd         = 0;
            for (int L = 0; L < DimLengths[0]; ++L) {
               for (int K = 0; K < DimLengths[1]; ++K) {
                  for (int J = 0; J < DimLengths[2]; ++J) {
                     for (int I = 0; I < DimLengths[3]; ++I) {
                        Data(L, K, J, I) = DataR8[VecAdd];
                        ++VecAdd;
                     }
                  }
               }
            }
            deepCopy(DataTmp, Data);
         }
         break;
      case 5:
         if (OnHost) {
            HostArray5DR8 Data = FieldPtr->getDataArray<HostArray5DR8>();
            int VecAdd         = 0;
            for (int M = 0; M < DimLengths[0]; ++M) {
               for (int L = 0; L < DimLengths[1]; ++L) {
                  for (int K = 0; K < DimLengths[2]; ++K) {
                     for (int J = 0; J < DimLengths[3]; ++J) {
                        for (int I = 0; I < DimLengths[4]; ++I) {
                           Data(M, L, K, J, I) = DataR8[VecAdd];
                           ++VecAdd;
                        }
                     }
                  }
               }
            }
         } else {
            Array5DR8 DataTmp  = FieldPtr->getDataArray<Array5DR8>();
            HostArray5DR8 Data = createHostMirrorCopy(DataTmp);
            int VecAdd         = 0;
            for (int M = 0; M < DimLengths[0]; ++M) {
               for (int L = 0; L < DimLengths[1]; ++L) {
                  for (int K = 0; K < DimLengths[2]; ++K) {
                     for (int J = 0; J < DimLengths[3]; ++J) {
                        for (int I = 0; I < DimLengths[4]; ++I) {
                           Data(M, L, K, J, I) = DataR8[VecAdd];
                           ++VecAdd;
                        }
                     }
                  }
               }
            }
            deepCopy(DataTmp, Data);
         }
         break;
      } // end switch NDims
      break; // end R8 fields

   default:
      ABORT_ERROR("IOStream readFieldData "
                  "Invalid data type while reading field {} for stream {}",
                  FieldName, Name);

   } // end switch data type

   // Clean up the decomp
   IO::destroyDecomp(DecompID);

   return Err;

} // End readFieldData

//------------------------------------------------------------------------------
// Reads a stream if it is time. This is the internal read function used by the
// public read interface.
Error IOStream::readStream(
    const Clock *ModelClock, // [in] model clock for getting time
    Metadata &ReqMetadata,   // [inout] global metadata to extract from file
    bool ForceRead           // [in] optional: read even if not time
) {

   Error Err; // returned error code

   // First check that this is an input stream
   if (Mode != IO::ModeRead)
      ABORT_ERROR("IOStream read: cannot read stream defined as output stream");

   // If it is not time to read, return
   if (!ForceRead) {
      if (!MyAlarm.isRinging() and !OnStartup)
         return Err;
      if (UseStartEnd) { // If time outside interval, return
         if (!StartAlarm.isRinging())
            return Err;
         if (EndAlarm.isRinging())
            return Err;
      }
   }

   // Get current simulation time and time string
   TimeInstant SimTime    = ModelClock->getCurrentTime();
   std::string SimTimeStr = SimTime.getString(5, 0, "_");

   // Reset alarms and flags
   if (OnStartup)
      OnStartup = false;
   if (MyAlarm.isRinging())
      MyAlarm.reset(SimTime);

   // Create filename
   std::string InFileName;
   // If using pointer files for this stream, read the filename from the pointer
   if (UsePointer) {
      std::ifstream PtrFile(PtrFilename);
      PtrFile >> InFileName;
      PtrFile.close();
   } else if (FilenameIsTemplate) {
      // create file name from template
      InFileName = buildFilename(Filename, SimTime);
   } else {
      InFileName = Filename;
   }

   // Open input file
   int InFileID;
   IO::openFile(InFileID, InFileName, Mode, IO::FmtDefault, ExistAction);

   // Read any requested global metadata
   for (auto Iter = ReqMetadata.begin(); Iter != ReqMetadata.end(); ++Iter) {
      std::string MetaName = Iter->first;
      std::any MetaTmp     = Iter->second;

      I4 MetaValI4;
      I8 MetaValI8;
      R4 MetaValR4;
      R8 MetaValR8;
      bool MetaValBool;
      std::string MetaValStr;

      // Read metadata based on type
      Error ErrRead;
      if (MetaTmp.type() == typeid(I8)) {
         ErrRead = IO::readMeta(MetaName, MetaValI8, InFileID, IO::GlobalID);
         ReqMetadata[MetaName] = MetaValI8;
      } else if (MetaTmp.type() == typeid(I4)) {
         ErrRead = IO::readMeta(MetaName, MetaValI4, InFileID, IO::GlobalID);
         ReqMetadata[MetaName] = MetaValI4;
      } else if (MetaTmp.type() == typeid(R8)) {
         ErrRead = IO::readMeta(MetaName, MetaValR8, InFileID, IO::GlobalID);
         ReqMetadata[MetaName] = MetaValR8;
      } else if (MetaTmp.type() == typeid(R4)) {
         ErrRead = IO::readMeta(MetaName, MetaValR4, InFileID, IO::GlobalID);
         ReqMetadata[MetaName] = MetaValR4;
      } else if (MetaTmp.type() == typeid(bool)) {
         // bool must be read as int
         ErrRead = IO::readMeta(MetaName, MetaValI4, InFileID, IO::GlobalID);
         if (MetaValI4 == 0) {
            MetaValBool = false;
         } else {
            MetaValBool = true;
         }
         ReqMetadata[MetaName] = MetaValBool;
      } else if (MetaTmp.type() == typeid(std::string)) {
         ErrRead = IO::readMeta(MetaName, MetaValStr, InFileID, IO::GlobalID);
         ReqMetadata[MetaName] = MetaValStr;
         // If ReqMetadata was initialized with a string literal, we detect
         // the type but replace it with a std::string
      } else if (MetaTmp.type() == typeid(const char *)) {
         ErrRead = IO::readMeta(MetaName, MetaValStr, InFileID, IO::GlobalID);
         ReqMetadata[MetaName] = MetaValStr;
      } else {
         ABORT_ERROR(
             "Metadata read failed: unknown data type for {} in file {}",
             MetaName, InFileName);
      }

      CHECK_ERROR_ABORT(ErrRead, "Error reading metadata {} from file {}",
                        MetaName, InFileName);

   } // end loop over requested metadata

   // Get dimensions from file and check that file has same dimension lengths
   std::map<std::string, int> AllDimIDs;
   defineAllDims(InFileID, AllDimIDs);

   // For each field in the contents, define field and read field data
   for (auto IFld = Contents.begin(); IFld != Contents.end(); ++IFld) {

      // Retrieve the field name and pointer
      std::string FieldName            = *IFld;
      std::shared_ptr<Field> ThisField = Field::get(FieldName);

      // Extract the data pointer and read the data array
      int FieldID; // not currently used but available if field metadata needed
      Err += readFieldData(ThisField, InFileID, AllDimIDs, FieldID);

   } // End loop over field list

   // Close input file
   IO::closeFile(InFileID);

   if (Err.isSuccess())
      LOG_INFO("Successfully read stream {} from file {}", Name, InFileName);

   // End of routine - return
   return Err;

} // End read

//------------------------------------------------------------------------------
// Writes stream. This is the internal member write function used by the
// public write interfaces.
void IOStream::writeStream(
    const Clock *ModelClock, // [in] Model clock needed for time stamps
    bool ForceWrite,         // [in] Optional: write even if not time
    bool FinalCall           // [in] Optional flag if called from finalize
) {

   Error Err; // default return code

   // First check that this is an output stream
   if (Mode != IO::ModeWrite)
      ABORT_ERROR("IOStream write: cannot write to an input stream");

   // If it is not time to write, return
   if (!ForceWrite) {
      bool StartupShutdown = OnStartup or (OnShutdown and FinalCall);
      if (!MyAlarm.isRinging() and !StartupShutdown)
         return;
      if (UseStartEnd) { // If time outside interval, return
         if (!StartAlarm.isRinging())
            return;
         if (EndAlarm.isRinging())
            return;
      }
   }

   // Get start time, current simulation time and time string
   TimeInstant StartTime  = ModelClock->getStartTime();
   TimeInstant SimTime    = ModelClock->getCurrentTime();
   std::string SimTimeStr = SimTime.getString(4, 0, "_");

   // Determine the time to use for the filename. The default is to
   // use the current time.
   TimeInstant FileTime = ModelClock->getCurrentTime();
   // For streams that need to write multiple frames or time slices per file,
   // and it is not time for a new file, then we use the time the file was
   // started which should be the last file alarm time.
   if (Multiframe and !FileAlarm.isRinging())
      FileTime = *(FileAlarm.getRingTimePrev());

   // Update the time field with elapsed time since simulation start
   TimeInterval ElapsedTime = SimTime - StartTime;
   R8 ElapsedTimeR8;
   ElapsedTime.get(ElapsedTimeR8, TimeUnits::Seconds);
   HostArray1DR8 OutTime("OutTime", 1);
   OutTime(0) = ElapsedTimeR8;
   Field::attachFieldData("time", OutTime);

   // Reset alarms and flags
   if (OnStartup)
      OnStartup = false;
   if (MyAlarm.isRinging())
      MyAlarm.reset(SimTime);
   if (FileAlarm.isRinging())
      FileAlarm.reset(SimTime);

   // Create filename
   std::string OutFileName;
   if (FilenameIsTemplate) {
      // create file name from template
      OutFileName = buildFilename(Filename, FileTime);
   } else {
      OutFileName = Filename;
   }

   // Open output file
   int OutFileID;
   IO::openFile(OutFileID, OutFileName, Mode, IO::FmtDefault, ExistAction);

   // For files with multiple frames or time slices, we need to determine the
   // default Frame number for time-dependent fields. If the frame/time already
   // exists in the file, we assume we are re-running the time period and
   // should over-write the existing frame. Otherwise, the new frame is the
   // next frame in the sequence.
   Frame = 0; // default is one frame in file
   if (Multiframe) {
      // Get the current number of frames in the file - the length of the
      // unlimited time dimension
      I4 TimeDimID;
      I4 NFrames;
      Err = IO::getDimFromFile(OutFileID, "time", TimeDimID, NFrames);
      if (Err.isFail() or NFrames == 0) {
         // If there is an error, we assume this is a new file in which
         // the dimension has not yet been written. Similarly, if NFrames is 0,
         // no frames have yet been written. In both cases, this is the first
         // frame in the file.
         Frame = 0;
      } else {
         // If there are frames in the file, we read the elapsed time for each
         // frame to determine whether the current slice already exists. If
         // equal (within roundoff), the slice exists and we overwrite with
         // the current frame. Otherwise, increment to the next frame
         R8 FrameTime;
         int TmpID;
         std::vector<int> TmpDimLengths; // empty dim length for scalar time
         for (int IFrame = 0; IFrame < NFrames; ++IFrame) {
            Err = IO::readNDVar(&FrameTime, "time", OutFileID, TmpID, IFrame,
                                &TmpDimLengths);
            CHECK_ERROR_ABORT(Err, "Error reading frame time in {}",
                              OutFileName);
            if (std::abs(ElapsedTimeR8 - FrameTime) < 1.e-5) { // overwrite
               Frame = IFrame;
               break;
            } else if (ElapsedTimeR8 > FrameTime) { // move to next frame
               Frame = IFrame + 1;
            } else { // time is earlier but doesn't match any frames
               ABORT_ERROR("Existing multiframe file {} appears to be using "
                           "different time intervals or is from the wrong "
                           "time period",
                           OutFileName);
            } // end if elapsed time matches
         } // end loop over existing frames
      } // end if nframes
   } // end if multiframe

   // Write Metadata for global metadata (Code and Simulation)
   // Only needs to be written for a new file
   if (Frame < 1) {
      writeFieldMeta(CodeMeta, OutFileID, IO::GlobalID);
      writeFieldMeta(SimMeta, OutFileID, IO::GlobalID);
   }

   // Create and write a field for any global data that is file or time
   // specific.
   std::shared_ptr<Field> FileField = Field::create("FileField");
   // Add the current simulation time
   std::string SimTimeName = "SimulationTime";
   FileField->addMetadata(SimTimeName, SimTimeStr);

   // If it is a multi-frame file, add the time for this frame
   if (Multiframe) {
      SimTimeName = SimTimeName + std::to_string(Frame);
      FileField->addMetadata(SimTimeName, SimTimeStr);
   }
   // Write and then destroy temporary field
   writeFieldMeta("FileField", OutFileID, IO::GlobalID);
   Field::destroy("FileField");

   // Assign dimension IDs for all defined dimensions
   std::map<std::string, int> AllDimIDs;
   defineAllDims(OutFileID, AllDimIDs);

   // Define each field and write field metadata
   std::map<std::string, int> FieldIDs;
   I4 NDims;
   std::vector<std::string> DimNames;
   std::vector<int> FieldDims;
   for (auto IFld = Contents.begin(); IFld != Contents.end(); ++IFld) {

      // Retrieve the field pointer
      std::string FieldName            = *IFld;
      std::shared_ptr<Field> ThisField = Field::get(FieldName);

      // Retrieve the dimensions for this field and determine dim IDs
      NDims = ThisField->getNumDims();
      if (NDims > 0) {
         DimNames.resize(NDims);
         FieldDims.resize(NDims);
         ThisField->getDimNames(DimNames);
      }
      // If this is a time-dependent field, we insert the unlimited time
      // dimension as the first dimension (for field definition only)
      if (ThisField->isTimeDependent()) {
         ++NDims;
         DimNames.insert(DimNames.begin(), "time");
         FieldDims.resize(NDims);
      }
      // Get the dim IDs
      for (int IDim = 0; IDim < NDims; ++IDim) {
         std::string DimName = DimNames[IDim];
         FieldDims[IDim]     = AllDimIDs[DimName];
      }

      // Determine the data type and convert to IODataType
      // Reduce floating point precision if requested
      IO::IODataType MyIOType = getFieldIOType(ThisField);

      // Define the field and assign a FieldID
      int FieldID =
          defineVar(OutFileID, FieldName, MyIOType, NDims, FieldDims.data());

      FieldIDs[FieldName] = FieldID;

      // Now we can write the field metadata
      if (Frame < 1) { // only write if it's the first time
         writeFieldMeta(FieldName, OutFileID, FieldID);
      }
   }

   // Now write data arrays for all fields in contents
   for (auto IFld = Contents.begin(); IFld != Contents.end(); ++IFld) {

      // Retrieve the field pointer and FieldID
      std::string FieldName            = *IFld;
      std::shared_ptr<Field> ThisField = Field::get(FieldName);
      int FieldID                      = FieldIDs[FieldName];

      // Extract and write the data array
      this->writeFieldData(ThisField, OutFileID, FieldID, AllDimIDs);
   }

   // Close output file
   IO::closeFile(OutFileID);

   // If using pointer files for this stream, write the filename to the pointer
   // after the file is successfully written
   if (UsePointer) {
      std::ofstream PtrFile(PtrFilename, std::ios::trunc);
      PtrFile << OutFileName << std::endl;
      PtrFile.close();
   }

   LOG_INFO("Successfully wrote stream {} to file {}", Name, OutFileName);

   // End of routine - return
   return;

} // end writeStream

//------------------------------------------------------------------------------
// Removes a single IOStream from the list of all streams.
void IOStream::erase(const std::string &StreamName // Name of IOStream to remove
) {
   AllStreams.erase(StreamName); // use the map erase function to remove
} // End erase

//------------------------------------------------------------------------------
// Private utility functions for read/write
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// Determines the IO Data type to use for a given field, taking into
// account the field's type and any reduced precision conversion
IO::IODataType IOStream::getFieldIOType(
    std::shared_ptr<Field> FieldPtr // [in] pointer to Field
) {

   IO::IODataType ReturnType;

   // Get the field data type
   ArrayDataType MyType = FieldPtr->getType();

   // Determine IO data type based on field type and any reduced precision
   // conversion
   switch (MyType) {
   case ArrayDataType::I4:
      ReturnType = IO::IOTypeI4;
      break;
   case ArrayDataType::I8:
      ReturnType = IO::IOTypeI8;
      break;
   case ArrayDataType::R4:
      ReturnType = IO::IOTypeR4;
      break;
   case ArrayDataType::R8:
      if (ReducePrecision and !FieldPtr->retainPrecision()) {
         ReturnType = IO::IOTypeR4;
      } else {
         ReturnType = IO::IOTypeR8;
      }
      break;
   default:
      std::string FieldName = FieldPtr->getName();
      ABORT_ERROR("IOStream getIOType: "
                  "Cannot determine data type for field {}",
                  FieldName);
      break;
   }

   return ReturnType;

} // end getIOType

//------------------------------------------------------------------------------
// Builds a filename based on time information and a filename template
// where special tokens are translated as:
//    $SimTime  = simulation time in form YYYY-MM-DD_hh.mm.ss
//    $WallTime = actual wallclock time in form YYYY-MM-DD_hh.mm.ss
// and individual components of simulation time can be used to customize
// the time stamp using:
//    $Y = year    part of simulation time stamp
//    $M = month   part of simulation time stamp
//    $D = day     part of simulation time stamp
//    $h = hour    part of simulation time stamp
//    $m = minute  part of simulation time stamp
//    $s = seconds part of simulation time stamp
std::string IOStream::buildFilename(
    const std::string &FilenameTemplate, // [in] template string for name
    const TimeInstant &FileTime          // [in] sim time to use for filename
) {

   // Start with input template
   std::string Outfile = FilenameTemplate;

   // Check if wallclock time is requested in the template - check for
   // multiple variations
   size_t Pos = Outfile.find("$WallTime");
   if (Pos == std::string::npos)
      Pos = Outfile.find("$Walltime");
   if (Pos == std::string::npos)
      Pos = Outfile.find("$walltime");
   if (Pos == std::string::npos)
      Pos = Outfile.find("$WALLTIME");
   // If wallclock time is requested, replace with wallclock string
   // in format YYYY-MM-DD_HH.MM.SS
   if (Pos != std::string::npos) { // wallclock string was found
      // Get wallclock time and convert to string
      std::time_t Walltime  = std::time(nullptr);
      std::tm *WalltimeInfo = std::localtime(&Walltime);
      char Walltimestamp[20];
      std::strftime(Walltimestamp, sizeof(Walltimestamp), "%Y-%m-%d_%H.%M.%S",
                    WalltimeInfo);
      std::string WalltimeStr = Walltimestamp; // convert to std::string
      // Replace template string with wallclock string
      Outfile.replace(Pos, 9, WalltimeStr);
   }

   // For many template options we also need the components of the current
   // sim time so we extract them here
   I8 IYear;
   I8 IMonth;
   I8 IDay;
   I8 IHour;
   I8 IMin;
   I8 ISecW;
   I8 ISecN;
   I8 ISecD;
   FileTime.get(IYear, IMonth, IDay, IHour, IMin, ISecW, ISecN, ISecD);

   // Convert these to strings (will pad with zeros if needed later)
   std::string SYear  = std::to_string(IYear);
   std::string SMonth = std::to_string(IMonth);
   std::string SDay   = std::to_string(IDay);
   std::string SHour  = std::to_string(IHour);
   std::string SMin   = std::to_string(IMin);
   std::string SSec   = std::to_string(ISecW);

   // Now check for the $SimTime token and replace with a sim time string
   Pos = Outfile.find("$SimTime");
   if (Pos == std::string::npos)
      Pos = Outfile.find("$Simtime");
   if (Pos == std::string::npos)
      Pos = Outfile.find("$simtime");
   if (Pos == std::string::npos)
      Pos = Outfile.find("$SIMTIME");
   if (Pos != std::string::npos) { // Found the sim time token
      // Convert the SimTime to a string
      int YearLength         = SYear.length();
      int MinWidth           = 4;
      int YearWidth          = std::max(YearLength, MinWidth);
      std::string SimTimeStr = FileTime.getString(YearWidth, 0, "_");
      Outfile.replace(Pos, 8, SimTimeStr);
   }

   // Check for each of the other standard tokens and replace with
   // appropriate strings with padding as necessary
   int NPad    = 0;
   int SLength = 0;

   // Year is requested. The length is assumed to be at least 4 and
   // padded as needed
   if ((Pos = Outfile.find("$Y")) != std::string::npos) {
      int SLength = SYear.length();
      if (SLength < 4) {
         NPad = 4 - SLength;
         SYear.insert(0, NPad, '0');
      }
      Outfile.replace(Pos, 2, SYear);
   }

   // Month is requested, pad to 2 digit month and insert into filename
   if ((Pos = Outfile.find("$M")) != std::string::npos) {
      SLength = SMonth.length();
      if (SLength < 2) {
         NPad = 2 - SLength;
         SMonth.insert(0, NPad, '0');
      }
      Outfile.replace(Pos, 2, SMonth);
   }

   // Day is requested, pad to 2 digit month and insert into filename
   if ((Pos = Outfile.find("$D")) != std::string::npos) {
      SLength = SDay.length();
      if (SLength < 2) {
         NPad = 2 - SLength;
         SDay.insert(0, NPad, '0');
      }
      Outfile.replace(Pos, 2, SDay);
   }

   // Hour is requested, pad to 2 digit hour and insert into filename
   if ((Pos = Outfile.find("$h")) != std::string::npos) {
      SLength = SHour.length();
      if (SLength < 2) {
         NPad = 2 - SLength;
         SHour.insert(0, NPad, '0');
      }
      Outfile.replace(Pos, 2, SHour);
   }

   // Minute is requested, pad to 2 digit minute and insert into filename
   if ((Pos = Outfile.find("$m")) != std::string::npos) {
      SLength = SMin.length();
      if (SLength < 2) {
         NPad = 2 - SLength;
         SMin.insert(0, NPad, '0');
      }
      Outfile.replace(Pos, 2, SMin);
   }

   // Second is requested. We have already rounded down to whole
   // seconds, so pad that to two digits if needed
   if ((Pos = Outfile.find("$s")) != std::string::npos) {
      SLength = SSec.length();
      if (SLength < 2) {
         NPad = 2 - SLength;
         SSec.insert(0, NPad, '0');
      }
      Outfile.replace(Pos, 2, SSec);
   }

   // Add netcdf suffix
   Outfile = Outfile + ".nc";

   return Outfile;

} // End buildFilename

//------------------------------------------------------------------------------
// Sets ReducePrecision flag based on an input string, performing string
// manipulation for case insensitive comparison

void IOStream::setPrecisionFlag(
    const std::string &PrecisionString ///< [in] precision from input YAML
) {

   // Set default value
   ReducePrecision = false;

   // Convert input string to lowercase for easier comparison
   std::string PrecComp = PrecisionString;
   std::transform(PrecComp.begin(), PrecComp.end(), PrecComp.begin(),
                  [](unsigned char c) { return std::tolower(c); });

   // Reduced precision if input is single
   if (PrecComp == "single") {
      ReducePrecision = true;
   } else if (PrecComp == "double") {
      ReducePrecision = false;
   } else {
      ABORT_ERROR("IOStream setPrecisionFlag: "
                  "Unknown precision {} in stream {}",
                  PrecisionString, Name);
   }

} // End setPrecisionFlag

} // namespace OMEGA
//===----------------------------------------------------------------------===//
