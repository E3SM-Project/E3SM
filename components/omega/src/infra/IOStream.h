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
/// Streams are added to the input configuration file as shown in this example
/// \ConfigInput
/// # IO Streams (examples for Default config)
/// IOStreams:
///    # InitialState should only be used when starting from scratch.
///    # For restart runs, the frequency units should be changed from
///    # "OnStartup" to "never" so that the initial state file is not read.
///    InitialState:
///      UsePointerFile: false
///      Filename: OmegaMesh.nc
///      Mode: read
///      Precision: double
///      Freq: 1
///      FreqUnits: OnStartup
///      UseStartEnd: false
///      Contents:
///        - State
///        - Base
///    # Restarts are used to initialize for all job submissions after the very
///    # first startup job. We use UseStartEnd with a start time just after the
///    # simulation start time so that omega does not attempt to use a restart
///    # for the first startup job.
///    RestartRead:
///      UsePointerFile: true
///      PointerFilename: ocn.pointer
///      Mode: read
///      Precision: double
///      Freq: 1
///      FreqUnits: OnStartup
///      UseStartEnd: true
///      StartTime: 0001-01-01_00:00:01
///      EndTime: 99999-12-31_00:00:00
///      Contents:
///        - Restart
///    # Sample restart output stream
///    RestartWrite:
///      UsePointerFile: true
///      PointerFilename: ocn.pointer
///      Filename: ocn.restart.$Y-$M-$D_$h.$m.$s
///      Mode: write
///      IfExists: replace
///      Precision: double
///      Freq: 6
///      FreqUnits: months
///      UseStartEnd: false
///      Contents:
///        - Restart
///    # Sample history file - values at the specified time
///    History:
///      UsePointerFile: false
///      Filename: ocn.hist.$SimTime
///      Mode: write
///      # File format is only needed if differs from default (pnetcdf)
///      # FileFormat: adios
///      IfExists: replace
///      Precision: double
///      Freq: 1
///      FreqUnits: months
///      UseStartEnd: false
///      Contents:
///        - Tracers
///        - State
///        - SshCellDefault
///    # Sample high-frequency output. Limited fields and time duration
///    Highfreq:
///      UsePointerFile: false
///      Filename: ocn.hifreq.$Y-$M
///      Mode: write
///      # File format is only needed if differs from default (pnetcdf)
///      # FileFormat: adios
///      IfExists: append
///      Precision: single
///      Freq: 10
///      FreqUnits: days
///      FileFreq: 1
///      FileFreqUnits: months
///      UseStartEnd: true
///      StartTime: 0001-06-01_00:00:00
///      EndTime: 0001-07-31_00:00:00
///      Contents:
///        - Tracers
/// \EndConfigInput
///
//===----------------------------------------------------------------------===//

#include "Config.h"
#include "DataTypes.h"
#include "Dimension.h"
#include "Error.h"
#include "Field.h"
#include "IO.h"
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
   IO::FileFmt FileFormat;   ///< file format (eg pnetcdf, hdf5, adios)

   IO::Mode Mode;        ///< mode (read or write)
   bool ReducePrecision; ///< flag to use 32-bit precision for 64-bit floats
   Alarm MyAlarm;        ///< time mgr alarm for read/write
   Alarm FileAlarm;      ///< time mgr alarm for opening a new file
   bool OnStartup;       ///< flag to read/write on model startup
   bool OnShutdown;      ///< flag to read/write on model shutdown
   bool Multiframe;      ///< flag for multiple frames/time slices in file
   int Frame;            ///< current frame/slice for multi-frame stream

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
   static void create(const std::string &StreamName, ///< [in] name of stream
                      Config &StreamConfig, ///< [in] stream configuration
                      Clock *&ModelClock    ///< [inout] Omega model clock
   );

   /// Read all dimensions from an input file and determine the dimension ID.
   /// The file must be in data mode.
   void readAllDims(
       int FileID, ///< [in] id assigned to the IO file
       std::map<std::string, int> &AllDimIDs ///< [out] dim name, assigned ID
   );

   /// Define all dimensions used. If some dimensions have already be read
   /// from the file (ie when appending to an existing file), only the dims
   /// not already read will be defined. The file must be in define mode.
   /// A map of all dimension IDs is returned.
   void defineAllDims(
       int FileID, ///< [in] id assigned to the IO file
       std::map<std::string, int> &AllDimIDs ///< [out] dim name, assigned ID
   );

   /// Computes the parallel decomposition (offsets) for a field.
   /// Needed for parallel I/O
   void computeDecomp(
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
   Error readStream(
       const Clock *ModelClock, ///< [in] Model clock for alarms, time stamp
       Metadata &ReqMetadata, ///< [inout] global metadata to extract from file
       bool ForceRead = false ///< [in] Optional: read even if not time
   );

   /// Private function that performs most of the stream write - called by the
   /// public write method
   void writeStream(
       const Clock *ModelClock, ///< [in] Model clock for alarms, time stamp
       bool ForceWrite = false, ///< [in] Optional: write even if not time
       bool FinalCall  = false  ///< [in] Optional flag for shutdown
   );

   /// Write all metadata associated with a field
   void writeFieldMeta(std::string FieldName, ///< [in] metadata from field;
                       int FileID,            ///< [in] id assigned to open file
                       int FieldID            ///< [in] id assigned to the field
   );

   /// Write a field's data array, performing any manipulations to reduce
   /// precision or move data between host and device
   void
   writeFieldData(std::shared_ptr<Field> FieldPtr, ///< [in] field to write
                  int FileID,  ///< [in] id assigned to open file
                  int FieldID, ///< [in] id assigned to the field
                  std::map<std::string, int> &AllDimIDs ///< [in] dimension IDs
   );

   /// Read a field's data array, performing any manipulations to reduce
   /// precision or move data between host and device
   Error
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
       const TimeInstant &FileTime          ///< [in] time to use in filename
   );

   /// Sets ReducePrecision flag based on an input string, performing string
   /// manipulation for case insensitive comparison
   void setPrecisionFlag(
       const std::string &PrecisionString ///< [in] precision from input YAML
   );

 public:
   //---------------------------------------------------------------------------
   /// Default empty constructor
   IOStream();

   //---------------------------------------------------------------------------
   /// Creates all streams defined in the input configuration file. This does
   /// not validate the contents of the streams since the relevant Fields
   /// may not have been defined yet.
   static void init(Clock *&ModelClock ///< [inout] Omega model clock
   );

   //---------------------------------------------------------------------------
   /// Overloaded init with no args, helpful for tests where no Clock exists
   static void init(void);

   //---------------------------------------------------------------------------
   /// Performs a final write of any streams that have the OnShutdown option and
   /// then removes all streams to clean up.
   static void
   finalize(const Clock *ModelClock ///< [in] Model clock needed for time stamps
   );

   //---------------------------------------------------------------------------
   /// Overloaded finalize with no args, helpful for test where no Clock exists
   static void finalize(void);

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
   /// Reads a stream if it is time.
   static Error read(const std::string &StreamName, ///< [in] Name of stream
                     const Clock *ModelClock, ///< [in] Model clock w time info
                     Metadata &ReqMetadata,   ///< [inout] Metadata desired
                     bool ForceRead = false ///< [in] opt: read even if not time
   );

   //---------------------------------------------------------------------------
   /// Overloaded simplified Read, reads stream regardless of time.
   static Error read(const std::string &StreamName); ///< [in] Name of stream

   //---------------------------------------------------------------------------
   /// Writes a stream if it is time.
   static void
   write(const std::string &StreamName, ///< [in] Name of stream
         const Clock *ModelClock,       ///< [in] Model clock for time stamps
         bool ForceWrite = false        ///< [in] opt: write even if not time
   );

   //---------------------------------------------------------------------------
   /// Loops through all streams and writes them if it is time. This is
   /// useful if most I/O is consolidated at one point (eg end of step).
   static void
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
