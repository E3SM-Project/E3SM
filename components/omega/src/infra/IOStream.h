#ifndef OMEGA_IOSTREAM_H
#define OMEGA_IOSTREAM_H
//===-- infra/IOStream.h - IO stream class ----------------------*- C++ -*-===//
//
/// \file
/// \brief Defines IO Stream class and methods
///
/// This header defines classes and methods for IO Streams. IOStreams define
/// reading and writing of fields from/to a data file. Each stream contains
/// information on the file location, the frequency of input/output, the
/// fields to be read/written and other information necessary for file I/O.
/// Note that this is different from the C++ stdlib iostreams
///
//===----------------------------------------------------------------------===//

#include "Config.h"
#include "DataTypes.h"
#include "Dimension.h"
#include "Field.h"
#include "IO.h"
#include "Logging.h"
#include "TimeMgr.h" // need Alarms, TimeInstant
#include <map>
#include <memory>
#include <set>
#include <string>

namespace OMEGA {

class IOStream {

 private:
   /// Store and maintain all defined streams
   static std::map<std::string, std::shared_ptr<IOStream>> AllStreams;

   /// Private variables specific to a stream
   std::string Name;         ///< name of stream
   std::string Filename;     ///< filename or filename template (with path)
   bool FilenameIsTemplate;  ///< true if the filename is a template
   IO::IfExists ExistAction; ///< action if file exists (write only)

   IO::Mode Mode;        ///< mode (read or write)
   bool ReducePrecision; ///< flag to use 32-bit precision for 64-bit floats
   Alarm MyAlarm;        ///< time mgr alarm for read/write
   bool OnStartup;       ///< flag to read/write on model startup
   bool OnShutdown;      ///< flag to read/write on model shutdown

   /// A pointer file is used if we wish OMEGA to read the name of the file
   /// from another file. This is useful for writing the name of a restart
   /// file to be picked up by the next job submitted so that the input
   /// configuration does not need to change with each submission
   bool UsePointer;         ///< flag for using a pointer file
   std::string PtrFilename; ///< name of pointer file

   /// Use a start and end time to define an interval in which stream is active
   /// The start is inclusive but the end time is not.
   bool UseStartEnd; ///< flag for using start, end times
   Alarm StartAlarm; ///< alarm ringing if equal to or past start time
   Alarm EndAlarm;   ///< alarm ringing if equal to or past end time

   /// Contents of stream in the form of a set of Field names
   std::set<std::string> Contents;

   /// Flag to determine whether the Contents have been validated or not
   bool Validated;

   //---- Private utility functions to support public interfaces
   /// Creates a new stream and adds to the list of all streams, based on
   /// options in the input model configuration. This routine is called by
   /// the IOStreams initialize function. It requires an initialized model
   /// clock so that stream alarm can be attached to this clock during creation.
   static int create(const std::string &StreamName, ///< [in] name of stream
                     Config &StreamConfig, ///< [in] input stream configuration
                     Clock *&ModelClock    ///< [inout] Omega model clock
   );

   /// Define all dimensions used. Returns an error code as well as a map
   /// of dimension names to defined dimension IDs.
   int defineAllDims(
       int FileID, ///< [in] id assigned to the IO file
       std::map<std::string, int> &AllDimIDs ///< [out] dim name, assigned ID
   );

   /// Computes the parallel decomposition (offsets) for a field.
   /// Needed for parallel I/O
   int computeDecomp(
       std::shared_ptr<Field> FieldPtr, ///< [in] field
       int &DecompID, ///< [out] ID assigned to the defined decomposition
       I4 &LocalSize, ///< [out] size of the local array for this field
       std::vector<int> &DimLengths // [out] local dim lengths
   );

   /// Retrieves field size information for non-distributed fields
   /// (distributed fields get this info from computeDecomp)
   void getFieldSize(std::shared_ptr<Field> FieldPtr, ///< [in] field
                     I4 &LocalSize, ///< [out] size of the local array
                     std::vector<int> &DimLengths ///< [out] local dim lengths
   );

   /// Private function that performs most of the stream read - called by the
   /// public read method
   int readStream(
       const Clock *ModelClock, ///< [in] Model clock for alarms, time stamp
       Metadata &ReqMetadata, ///< [inout] global metadata to extract from file
       bool ForceRead = false ///< [in] Optional: read even if not time
   );

   /// Private function that performs most of the stream write - called by the
   /// public write method
   int writeStream(
       const Clock *ModelClock, ///< [in] Model clock for alarms, time stamp
       bool ForceWrite = false, ///< [in] Optional: write even if not time
       bool FinalCall  = false  ///< [in] Optional flag for shutdown
   );

   /// Write all metadata associated with a field
   int writeFieldMeta(std::string FieldName, ///< [in] metadata from field;
                      int FileID,            ///< [in] id assigned to open file
                      int FieldID            ///< [in] id assigned to the field
   );

   /// Write a field's data array, performing any manipulations to reduce
   /// precision or move data between host and device
   int
   writeFieldData(std::shared_ptr<Field> FieldPtr, ///< [in] field to write
                  int FileID,  ///< [in] id assigned to open file
                  int FieldID, ///< [in] id assigned to the field
                  std::map<std::string, int> &AllDimIDs ///< [in] dimension IDs
   );

   /// Read a field's data array, performing any manipulations to reduce
   /// precision or move data between host and device
   int
   readFieldData(std::shared_ptr<Field> FieldPtr, ///< [in] field to read
                 int FileID, ///< [in] id assigned to open file
                 std::map<std::string, int> &AllDimIDs, ///< [in] dimension IDs
                 int &FieldID ///< [out] id assigned to the field
   );

   /// Determines the IO Data type to use for a given field, taking into
   /// account the field's type and any reduced precision conversion
   IO::IODataType getFieldIOType(
       std::shared_ptr<Field> FieldPtr ///< [in] field to extract type
   );

   /// Builds a filename based on time information and a filename template
   /// where special tokens are translated as:
   ///    $SimTime  = simulation time in form YYYY-MM-DD_hh.mm.ss
   ///    $WallTime = actual wallclock time in form YYYY-MM-DD_hh.mm.ss
   /// and individual components of simulation time can be used to customize
   /// the time stamp using:
   ///    $Y = year    part of simulation time stamp
   ///    $M = month   part of simulation time stamp
   ///    $D = day     part of simulation time stamp
   ///    $h = hour    part of simulation time stamp
   ///    $m = minute  part of simulation time stamp
   ///    $s = seconds part of simulation time stamp
   static std::string buildFilename(
       const std::string &FilenameTemplate, ///< [in] template string for name
       const Clock *ModelClock              ///< [in] model clock for sim time
   );

   /// Sets ReducePrecision flag based on an input string, performing string
   /// manipulation for case insensitive comparison
   void setPrecisionFlag(
       const std::string &PrecisionString ///< [in] precision from input YAML
   );

 public:
   //---------------------------------------------------------------------------
   /// Return codes - these will be removed once Error Handler is completed

   static constexpr int Success{0}; ///< Successful read/write completion
   static constexpr int Skipped{1}; ///< Normal early return (eg if not time)
   static constexpr int Fail{2};    ///< Fail

   //---------------------------------------------------------------------------
   /// Default empty constructor
   IOStream();

   //---------------------------------------------------------------------------
   /// Creates all streams defined in the input configuration file. This does
   /// not validate the contents of the streams since the relevant Fields
   /// may not have been defined yet. Returns an error code.
   static int init(Clock *&ModelClock ///< [inout] Omega model clock
   );

   //---------------------------------------------------------------------------
   /// Performs a final write of any streams that have the OnShutdown option and
   /// then removes all streams to clean up. Returns an error code.
   static int
   finalize(const Clock *ModelClock ///< [in] Model clock needed for time stamps
   );

   //---------------------------------------------------------------------------
   /// Retrieves a previously defined stream by name.
   static std::shared_ptr<IOStream>
   get(const std::string &StreamName ///< [in] name of stream to retrieve
   );

   //---------------------------------------------------------------------------
   /// Adds a field to the contents of a stream. Because streams may be created
   /// before all Fields have been defined, we only store the name. Validity
   /// is either checked during read/write or can be checked using the validate
   /// function.
   void addField(const std::string &FieldName ///< [in] Name of field
   );

   //---------------------------------------------------------------------------
   /// Removes a field from the contents. Provided for symmetry, but not
   /// often used.
   void removeField(const std::string &FieldName ///< [in] Name of field
   );

   //---------------------------------------------------------------------------
   /// Validates an IOStream and its contents. If used, this must be called at
   /// the end of initialization to ensure all Fields have been defined.
   /// This function also replaces any group names by the list of field members
   /// so that they can individually be checked and later read/write functions
   /// have the complete list of fields.  Returns true if all contents and
   /// variables are valid.
   bool validate();

   //---------------------------------------------------------------------------
   /// Determines whether the contents of the stream have been validated
   /// or not.
   bool isValidated();

   //---------------------------------------------------------------------------
   /// Validates all streams and their contents. If used, this must be called at
   /// the end of initialization to ensure all Fields have been defined.
   /// Returns true if all streams are valid.
   static bool validateAll();

   //---------------------------------------------------------------------------
   /// Reads a stream if it is time. Returns an error code.
   static int read(const std::string &StreamName, ///< [in] Name of stream
                   const Clock *ModelClock, ///< [in] Model clock for time info
                   Metadata &ReqMetadata, ///< [inout] Metadata desired in file
                   bool ForceRead = false ///< [in] opt: read even if not time
   );

   //---------------------------------------------------------------------------
   /// Writes a stream if it is time. Returns an error code.
   static int
   write(const std::string &StreamName, ///< [in] Name of stream
         const Clock *ModelClock,       ///< [in] Model clock for time stamps
         bool ForceWrite = false        ///< [in] opt: write even if not time
   );

   //---------------------------------------------------------------------------
   /// Loops through all streams and writes them if it is time. This is
   /// useful if most I/O is consolidated at one point (eg end of step).
   static int
   writeAll(const Clock *ModelClock ///< [in] Model clock for time stamps
   );

   //---------------------------------------------------------------------------
   /// Removes a single IOStream from the list of all streams.
   /// That process also decrements the reference counters for the
   /// shared pointers and removes them if those counters reach 0.
   static void erase(const std::string &StreamName /// Name of stream to remove
   );

}; // end class IOStream

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
#endif // OMEGA_IOSTREAM_H
