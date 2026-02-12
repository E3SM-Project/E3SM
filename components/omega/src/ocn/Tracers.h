#ifndef OMEGA_TRACERS_H
#define OMEGA_TRACERS_H
//===-- Tracers.h - tracer class definition ---------------------*- C++ -*-===//
//
/// \file
/// \brief Define Tracers class
///
/// The Tracers class define tracer arrays for host and device memories and
/// methods to handle the arrays. All of variables and methods are static
/// because once tracers are created, they never change.
/// The full list of tracers is defined in an include file TracerDefs.inc
/// and the specific tracers to use in a simulation are selected from that
/// list using the input configuration:
/// \ConfigInput
/// # List of requested tracers for a run
///  Tracers:
///    # Tracer group, list of individual tracers
///    Base: [Temperature, Salinity]
///    Debug: [Debug1, Debug2, Debug3]
/// \EndConfigInput
//
//===----------------------------------------------------------------------===//

#include "DataTypes.h"
#include "Decomp.h"
#include "Field.h"
#include "Halo.h"
#include "VertCoord.h"

namespace OMEGA {

/// The Tracers class provides a container for tracer arrays and methods
class Tracers {

 private:
   static I4 NumTracers;  ///< Total number of tracers defined at intialization
   static I4 NTimeLevels; ///< Number of time levels in tracer variable arrays
   static I4 NVertLayers; ///< Number of vertical layers in tracer arrays
   static I4 NCellsOwned; ///< Number of cells owned by this task
   static I4 NCellsAll;   ///< Total number of local cells (owned + all halo)
   static I4 NCellsSize;  ///< Array size (incl padding, bndy cell)

   inline static I4 IndxInvalid = -1;

   /// static storage of the tracer arrays for device and host memory spaces
   /// The arrays are 3 dimensional arrays with Tracer, Cell, and Vertical
   /// dimensions in order
   static std::vector<Array3DReal>
       TracerArrays; ///< TimeLevels -> [Tracer, Cell, Vert]
   static std::vector<HostArray3DReal>
       TracerArraysH; ///< TimeLevels -> [Tracer, Cell, Vert]

   /// maps for managing tracer groups
   /// Key of this map is a group name and
   /// Value is a pair of GroupStartIndex and GroupLength
   static std::map<std::string, std::pair<I4, I4>> TracerGroups;

   /// maps for matching tracer names with indices (both directions)
   static std::map<std::string, I4> TracerIndexes;
   static std::map<I4, std::string> TracerNames;

   /// for halo exchange
   static Halo *MeshHalo;

   /// Tracer dimension names
   static std::vector<std::string> TracerDimNames;

   /// Current time index
   /// this index is circular so that it returns to index 0
   /// if it is over max index
   static I4 CurTimeIndex; ///< Time dimension array index for current level

   /// get the time level index
   static I4 getTimeIndex(const I4 TimeLevel);

   /// locally defines all tracers but do not allocates memory
   static I4
   define(const std::string &Name,        ///< [in] Name of tracer
          const std::string &Description, ///< [in] Long name or description
          const std::string &Units,       ///< [in] Units
          const std::string &StdName,     ///< [in] CF standard Name
          const Real ValidMin,            ///< [in] min valid field value
          const Real ValidMax,            ///< [in] max valid field value
          const Real FillValue,           ///< [in] value for undef entries
          I4 &Index = IndxInvalid         ///< [out] (optional) index value
   );

 public:
   //---------------------------------------------------------------------------
   // Initialization
   //---------------------------------------------------------------------------
   /// read tracer defintions, allocate tracer arrays and initializes the
   /// tracers
   static void init();

#include "TracerDefs.inc"

   /// deallocates tracer arrays
   static I4 clear();

   //---------------------------------------------------------------------------
   // Query tracers
   //---------------------------------------------------------------------------

   /// get total number of tracers
   static I4 getNumTracers();

   /// get tracer index from tracer name
   static I4 getIndex(I4 &TracerIndex,              ///< [out] tracer index
                      const std::string &TracerName ///< [in] tracer name
   );

   /// get tracer name from tracer index
   static I4 getName(std::string &TracerName, ///< [out] tracer name
                     const I4 TracerIndex     ///< [in] tracer index
   );

   /// get a device array for all tracers
   static Array3DReal getAll(const I4 TimeLevel ///< [in] time level index
   );

   /// get a device array by tracer index
   static Array2DReal
   getByIndex(const I4 TimeLevel,  ///< [in] time level index
              const I4 TracerIndex ///< [in] global tracer index
   );

   /// get a device array by tracer name
   static Array2DReal
   getByName(const I4 TimeLevel,           ///< [in] time level index
             const std::string &TracerName ///< [in] global tracer name
   );

   /// get a host array for all tracers
   static HostArray3DReal
   getAllHost(const I4 TimeLevel ///< [in] time level index
   );

   /// get a host array by tracer index
   static HostArray2DReal
   getHostByIndex(const I4 TimeLevel,  ///< [in] time level index
                  const I4 TracerIndex ///< [in] global tracer index
   );

   /// get a host array by tracer name
   static HostArray2DReal getHostByName( ///< [out] tracer host array
       const I4 TimeLevel,               ///< [in] time level index
       const std::string &TracerName     ///< [in] global tracer name
   );

   /// get a field by tracer index. If not found, returns nullptr
   static std::shared_ptr<Field>
   getFieldByIndex(const I4 TracerIndex ///< [in] global tracer index
   );

   /// get a field by tracer name. If not found, returns nullptr
   static std::shared_ptr<Field>
   getFieldByName(const std::string &TracerName ///< [in] global tracer name
   );

   //---------------------------------------------------------------------------
   // Tracer group query
   //---------------------------------------------------------------------------

   /// return a vector of group names.
   static std::vector<std::string> getGroupNames();

   /// get a pair of (group start index, group length)
   static I4 getGroupRange(std::pair<I4, I4> &GroupRange, ///< [out] group range
                           const std::string &GroupName   ///< [in] group name
   );

   /// check if a tracer is a member of group by tracer index
   static bool
   isGroupMemberByIndex(const I4 TracerIndex,       ///< [in] tracer index
                        const std::string GroupName ///< [in] group name
   );

   /// check if a tracer is a member of group by tracer name
   static bool isGroupMemberByName(
       const std::string &TracerName, ///< [in] global tracer name
       const std::string &GroupName   ///< [in] group name
   );

   //---------------------------------------------------------------------------
   // Halo exchange and update time level
   //---------------------------------------------------------------------------

   /// Exchange halo
   static I4 exchangeHalo(const I4 TimeLevel ///< [in] tracer time level
   );

   /// increment time levels
   static void updateTimeLevels();

   //---------------------------------------------------------------------------
   // Device-Host data movement
   //---------------------------------------------------------------------------

   /// Copy tracers variables from host to device
   static void copyToDevice(const I4 TimeLevel ///< [in] tracer time level
   );

   /// Copy tracers variables from device to host
   static void copyToHost(const I4 TimeLevel ///< [in] tracer time level
   );

   //---------------------------------------------------------------------------
   // Forbid copy and move construction
   //---------------------------------------------------------------------------

   Tracers(const Tracers &) = delete;
   Tracers(Tracers &&)      = delete;

}; // end class Tracers

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
#endif // defined OMEGA_TRACERS_H
