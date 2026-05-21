#ifndef OMEGA_OCEANSTATE_H
#define OMEGA_OCEANSTATE_H
//===-- ocn/OceanState.h - ocean state --------------------*- C++ -*-===//
//
/// \file
/// \brief Contains the state variables for an OMEGA sub-domain
///
/// The OceanState class contains the non-tracer prognostic variable data for a
/// sub-domain of the global horizontal mesh.
//
//===----------------------------------------------------------------------===//

#include "DataTypes.h"
#include "Decomp.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "MachEnv.h"

#include <string>

namespace OMEGA {

/// A class for the ocean prognostic variable information

/// The OceanState class provides a container for the pseudo-thickness,
/// and normal velocity variables. It contains methods which handle
/// IO and time level updates.
class OceanState {

 private:
   void defineFields();

   Halo *MeshHalo;

   static OceanState *DefaultOceanState;

   static std::map<std::string, std::unique_ptr<OceanState>> AllOceanStates;

   /// Construct a new local state for a given decomposition
   OceanState(const std::string &Name, ///< [in] Name for mesh
              HorzMesh *Mesh,          ///< [in] Horizontal mesh
              Halo *MeshHalo_,         ///< [in] Halo for Mesh
              const int NVertLayers_,  ///< [in] Number of vertical layers
              const int NTimeLevels_   ///< [in] Number of time levels
   );

   // Forbid copy and move construction
   OceanState(const OceanState &) = delete;
   OceanState(OceanState &&)      = delete;

   // Current time index
   // this index is circular so that it returns to index 0
   // if it is over max index
   I4 CurTimeIndex; ///< Time dimension array index for current level

   /// Get the current time level index associated with a time level
   I4 getTimeIndex(const I4 TimeLevel) const;

 public:
   // Variables
   // Since these are used frequently, we make them public to reduce the
   // number of retrievals required.

   std::string Name;

   // Sizes and global IDs
   // Note that all sizes are actual counts (1-based) so that loop extents
   // should always use the 0:NCellsXX-1 form.

   I4 NCellsOwned; ///< Number of cells owned by this task
   I4 NCellsAll;   ///< Total number of local cells (owned + all halo)
   I4 NCellsSize;  ///< Array size (incl padding, bndy cell) for cell arrays

   I4 NEdgesOwned; ///< Number of edges owned by this task
   I4 NEdgesAll;   ///< Total number (owned+halo) of local edges
   I4 NEdgesSize;  ///< Array length (incl padding, bndy) for edge dim

   static const I4 MaxTimeLevels = 5; ///< Maximum number of time levels

   I4 NTimeLevels; ///< Number of time levels in state variable arrays
   I4 NVertLayers; ///< Number of vertical layers in state variable arrays

   // Prognostic variables

   std::vector<Array2DReal> PseudoThickness; ///< Device PseudoThickness array
   std::vector<HostArray2DReal>
       PseudoThicknessH; ///< Host PseudoThickness array

   std::vector<Array2DReal> NormalVelocity; ///< Device NormalVelocity array
   std::vector<HostArray2DReal> NormalVelocityH; ///< Host NormalVelocity array

   // Field names
   // These are appended with the State name for non-Default state instances
   std::string PseudoThicknessFldName; ///< Field name for PseudoThickness
   std::string NormalVelocityFldName;  ///< Field name for NormalVelocity
   std::string StateGroupName;

   // Methods

   /// Initialize Omega local state
   static int init();

   /// Create a new state by calling the constructor and put it in the
   /// AllOceanStates map
   static OceanState *
   create(const std::string &Name, ///< [in] Name for mesh
          HorzMesh *Mesh,          ///< [in] Horizontal mesh
          Halo *MeshHalo,          ///< [in] Halo for Mesh
          const int NVertLayers,   ///< [in] Number of vertical layers
          const int NTimeLevels    ///< [in] Number of time levels
   );

   /// Get pseudo-thickness device array at given time level
   Array2DReal getPseudoThickness(const I4 TimeLevel) const;

   /// Get pseudo-thickness host array at given time level
   HostArray2DReal getPseudoThicknessH(const I4 TimeLevel) const;

   /// Get normal velocity device array at given time level
   Array2DReal getNormalVelocity(const I4 TimeLevel) const;

   /// Get normal velocity host array at given time level
   HostArray2DReal getNormalVelocityH(const I4 TimeLevel) const;

   /// Exchange halo
   void exchangeHalo(const I4 TimeLevel);

   /// Swap time levels to update state arrays
   void updateTimeLevels();

   /// Copy state variables from host to device
   void copyToDevice(const I4 TimeLevel);

   /// Copy state variables from device to host
   void copyToHost(const I4 TimeLevel);

   /// Destructor - deallocates all memory and deletes an OceanState
   ~OceanState();

   /// Deallocates arrays
   static void clear();

   /// Remove state by name
   static void erase(std::string InName ///< [in] name of state to remove
   );

   static OceanState *getDefault();

   static OceanState *get(std::string name);

}; // end class OceanState

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
#endif // defined OMEGA_OCEANSTATE_H
