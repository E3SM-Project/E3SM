#ifndef OMEGA_TRACERS_H
#define OMEGA_TRACERS_H
//===-- ocn/Tracers.h - tracers --------------------*- C++ -*-===//
//
/// \file
/// \brief Define Tracers class
///
/// Tracers class define tracer arrays for host and device memories and
/// methods to handle the arrays. All of variables and methods are static
/// because once tracers are created, they never change.
//===----------------------------------------------------------------------===//

#include "DataTypes.h"
#include "Decomp.h"
#include "Field.h"
#include "Halo.h"

namespace OMEGA {

/// The Tracers class provides a container for tracer arrays and methods
class Tracers {

 private:
   // total number of tracers defined at intialization
   static I4 NumTracers;

   // static storage of the tracer arrays for device and host memory spaces
   // The arrays are 3 dimensional arrays with Tracer, Cell, and Vertical
   // dimensions in order
   static std::vector<Array3DReal>
       TracerArrays; ///< TimeLevels -> [Tracer, Cell, Vert]
   static std::vector<HostArray3DReal>
       TracerArraysH; ///< TimeLevels -> [Tracer, Cell, Vert]

   // maps for managing tracer groups
   // Key of this map is a group name and
   // Value is a pair of GroupStartIndex and GroupLength
   static std::map<std::string, std::pair<I4, I4>> TracerGroups;

   // maps for matching tracer names with indices (both directions)
   static std::map<std::string, I4> TracerIndexes;
   static std::map<I4, std::string> TracerNames;

   // for halo exchange
   static Halo *MeshHalo;

   // Tracer dimension names
   static std::vector<std::string> TracerDimNames;

   // Current time index
   // this index is circular so that it returns to index 0
   // if it is over max index
   static I4 CurTimeIndex; ///< Time dimension array index for current level

   // locally defines all tracers but do not allocates memory
   static I4
   define(const std::string &Name,        ///< [in] Name of tracer
          const std::string &Description, ///< [in] Long name or description
          const std::string &Units,       ///< [in] Units
          const std::string &StdName,     ///< [in] CF standard Name
          const Real ValidMin,            ///< [in] min valid field value
          const Real ValidMax,            ///< [in] max valid field value
          const Real FillValue            ///< [in] value for undef entries
   );

 public:
   static I4 NTimeLevels; ///< Number of time levels in tracer variable arrays
   static I4
       NVertLevels; ///< Number of vertical levels in tracer variable arrays
   static I4 NCellsOwned; ///< Number of cells owned by this task
   static I4 NCellsAll;   ///< Total number of local cells (owned + all halo)
   static I4
       NCellsSize; ///< Array size (incl padding, bndy cell) for cell arrays

   //---------------------------------------------------------------------------
   // Initialization
   //---------------------------------------------------------------------------
   /// read tracer defintions, allocate tracer arrays and initializes the
   /// tracers
   static I4 init();

   /// deallocates tracer arrays
   static I4 clear();

   //---------------------------------------------------------------------------
   // Query tracers
   //---------------------------------------------------------------------------

   // get total number of tracers
   static I4 getNumTracers();

   // get tracer index from tracer name
   static I4 getIndex(I4 &TracerIndex,              ///< [out] tracer index
                      const std::string &TracerName ///< [in] tracer name
   );

   // get tracer name from tracer index
   static I4 getName(std::string &TracerName, ///< [out] tracer name
                     const I4 TracerIndex     ///< [in] tracer index
   );

   // get a device array for all tracers
   static I4 getAll(Array3DReal &TracerArray, ///< [out] tracer device array
                    const I4 TimeLevel        ///< [in] time level index
   );

   // get a device array by tracer index
   static I4 getByIndex(Array2DReal &TracerArray, ///< [out] tracer device array
                        const I4 TimeLevel,       ///< [in] time level index
                        const I4 TracerIndex      ///< [in] global tracer index
   );

   // get a device array by tracer name
   static I4
   getByName(Array2DReal &TracerArray,     ///< [out] tracer device array
             const I4 TimeLevel,           ///< [in] time level index
             const std::string &TracerName ///< [in] global tracer name
   );

   // get a host array for all tracers
   static I4
   getAllHost(HostArray3DReal &TracerArrayH, ///< [out] tracer host array
              const I4 TimeLevel             ///< [in] time level index
   );

   // get a host array by tracer index
   static I4
   getHostByIndex(HostArray2DReal &TracerArrayH, ///< [out] tracer host array
                  const I4 TimeLevel,            ///< [in] time level index
                  const I4 TracerIndex           ///< [in] global tracer index
   );

   // get a host array by tracer name
   static I4
   getHostByName(HostArray2DReal &TracerArrayH, ///< [out] tracer host array
                 const I4 TimeLevel,            ///< [in] time level index
                 const std::string &TracerName  ///< [in] global tracer name
   );

   // get a field by tracer index. If not found, returns nullptr
   static std::shared_ptr<Field>
   getFieldByIndex(const I4 TracerIndex ///< [in] global tracer index
   );

   // get a field by tracer name. If not found, returns nullptr
   static std::shared_ptr<Field>
   getFieldByName(const std::string &TracerName ///< [in] global tracer name
   );

   //---------------------------------------------------------------------------
   // Tracer group query
   //---------------------------------------------------------------------------

   // return a vector of group names.
   static std::vector<std::string> getGroupNames();

   // get a pair of (group start index, group length)
   static I4 getGroupRange(std::pair<I4, I4> &GroupRange, ///< [out] group range
                           const std::string &GroupName   ///< [in] group name
   );

   // check if a tracer is a member of group by tracer index
   static bool
   isGroupMemberByIndex(const I4 TracerIndex,       ///< [in] tracer index
                        const std::string GroupName ///< [in] group name
   );

   // check if a tracer is a member of group by tracer name
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
   static I4 updateTimeLevels();

   //---------------------------------------------------------------------------
   // Device-Host data movement
   //---------------------------------------------------------------------------

   /// Copy tracers variables from host to device
   static I4 copyToDevice(const I4 TimeLevel ///< [in] tracer time level
   );

   /// Copy tracers variables from device to host
   static I4 copyToHost(const I4 TimeLevel ///< [in] tracer time level
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
